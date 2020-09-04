#include "s3.h"
using namespace himan;

#ifdef HAVE_S3
#include "debug.h"
#include "timer.h"
#include "util.h"
#include <iostream>
#include <libs3.h>
#include <mutex>
#include <string.h>  // memcpy

static std::once_flag oflag;

const char* access_key = 0;
const char* secret_key = 0;
const char* security_token = 0;

thread_local S3Status statusG = S3StatusOK;

void CheckS3Error(S3Status errarg, const char* file, const int line);

#define S3_CHECK(errarg) CheckS3Error(errarg, __FILE__, __LINE__)

std::vector<std::string> GetBucketAndFileName(const std::string& fullFileName)
{
	std::vector<std::string> ret;

	auto fileName = fullFileName;

	// strip protocol from string if it's there
	const auto pos = fullFileName.find("s3://");

	if (pos != std::string::npos)
	{
		fileName = fileName.erase(pos, 5);
	}

	// erase forward slash if exists (s3 buckets can't start with /)
	if (fileName[0] == '/')
	{
		fileName = fileName.erase(0, 1);
	}

	auto tokens = util::Split(fileName, "/");

	ret.push_back(tokens[0]);
	tokens.erase(std::begin(tokens), std::begin(tokens) + 1);

	std::string key;
	for (const auto& piece : tokens)
	{
		if (!key.empty())
		{
			key += "/";
		}
		key += piece;
	}

	ret.push_back(key);
	return ret;
}

inline void CheckS3Error(S3Status errarg, const char* file, const int line)
{
	if (errarg)
	{
		std::cerr << "Error at " << file << "(" << line << "): " << S3_get_status_name(errarg) << std::endl;
		himan::Abort();
	}
}

S3Status responsePropertiesCallback(const S3ResponseProperties* properties, void* callbackData)
{
	return S3StatusOK;
}

static void responseCompleteCallback(S3Status status, const S3ErrorDetails* error, void* callbackData)
{
	statusG = status;
	return;
}

thread_local S3ResponseHandler responseHandler = {&responsePropertiesCallback, &responseCompleteCallback};

static S3Status getObjectDataCallback(int bufferSize, const char* buffer, void* callbackData)
{
	himan::buffer* ret = static_cast<himan::buffer*>(callbackData);

	ret->data = static_cast<unsigned char*>(realloc(ret->data, ret->length + bufferSize));
	memcpy(ret->data + ret->length, buffer, bufferSize);
	ret->length += bufferSize;

	return S3StatusOK;
}

void Initialize()
{
	call_once(oflag, [&]() {
		access_key = getenv("S3_ACCESS_KEY_ID");
		secret_key = getenv("S3_SECRET_ACCESS_KEY");
		security_token = getenv("S3_SESSION_TOKEN");

		logger logr("s3");

		if (!access_key)
		{
			logr.Fatal("Environment variable S3_ACCESS_KEY_ID not defined");
			himan::Abort();
		}

		if (!secret_key)
		{
			logr.Fatal("Environment variable S3_SECRET_ACCESS_KEY not defined");
			himan::Abort();
		}

		S3_CHECK(S3_initialize("s3", S3_INIT_ALL, NULL));
	});
}

buffer s3::ReadFile(const file_information& fileInformation)
{
	Initialize();
	logger logr("s3");

	S3GetObjectHandler getObjectHandler = {responseHandler, &getObjectDataCallback};

	const auto bucketAndFileName = GetBucketAndFileName(fileInformation.file_location);
	const auto bucket = bucketAndFileName[0];
	const auto key = bucketAndFileName[1];

	buffer ret;

#ifdef S3_DEFAULT_REGION

	std::string region;

	if (fileInformation.file_server.find("amazonaws.com") != std::string::npos)
	{
		// extract region name from host name, assuming aws hostname like
		// s3.us-east-1.amazonaws.com

		auto tokens = util::Split(fileInformation.file_server, ".");

		if (tokens.size() != 4)
		{
			logr.Fatal("Hostname does not follow pattern s3.<regionname>.amazonaws.com");
		}
		else
		{
			logr.Trace("s3 authentication hostname: " + tokens[1]);
			region = tokens[1];
		}
	}

	const char* regionp = (region.empty()) ? nullptr : region.c_str();

	// clang-format off

        S3BucketContext bucketContext =
        {
                fileInformation.file_server.c_str(),
                bucket.c_str(),
                S3ProtocolHTTP,
                S3UriStylePath,
                access_key,
                secret_key,
                security_token,
                regionp
        };
#else

	// clang-format off

	S3BucketContext bucketContext = 
	{
		fileInformation.file_server.c_str(),
		bucket.c_str(),
		S3ProtocolHTTP,
		S3UriStylePath,
		access_key,
		secret_key,
		security_token
	};

	// clang-format on
#endif

	int count = 0;
	do
	{
		if (count > 0)
		{
			sleep(2 * count);
		}
		const unsigned long offset = fileInformation.offset.get();
		const unsigned long length = fileInformation.length.get();

		logr.Debug("Reading from host=" + fileInformation.file_server + " bucket=" + bucket + " key=" + key + " " +
		           std::to_string(offset) + ":" + std::to_string(length) + " (" + S3_get_status_name(statusG) + ")");

#ifdef S3_DEFAULT_REGION
		S3_get_object(&bucketContext, key.c_str(), NULL, offset, length, NULL, 0, &getObjectHandler, &ret);
#else
		S3_get_object(&bucketContext, key.c_str(), NULL, offset, length, NULL, &getObjectHandler, &ret);
#endif
		count++;
	} while (S3_status_is_retryable(statusG) && count < 3);

	switch (statusG)
	{
		case S3StatusOK:
			break;
		case S3StatusInternalError:
			logr.Error(std::string(S3_get_status_name(statusG)) + ": is there a proxy blocking the connection?");
			throw himan::kFileDataNotFound;
		case S3StatusFailedToConnect:
			logr.Error(std::string(S3_get_status_name(statusG)) + ": is proxy required but not set?");
			throw himan::kFileDataNotFound;
		case S3StatusErrorInvalidAccessKeyId:
			logr.Error(std::string(S3_get_status_name(statusG)) +
			           ": are Temporary Security Credentials used without security token (env: S3_SESSION_TOKEN)?");
			throw himan::kFileDataNotFound;
		default:
			logr.Error(S3_get_status_name(statusG));
			throw himan::kFileDataNotFound;
	}

	if (ret.length == 0)
	{
		throw himan::kFileDataNotFound;
	}

	return ret;
}

struct write_data
{
	himan::buffer buffer;
	size_t write_ptr;
};

static int putObjectDataCallback(int bufferSize, char* buffer, void* callbackData)
{
	write_data* data = static_cast<write_data*>(callbackData);
	int bytesWritten = 0;

	if (data->buffer.length)
	{
		bytesWritten =
		    static_cast<int>((static_cast<int>(data->buffer.length) > bufferSize) ? bufferSize : data->buffer.length);
		memcpy(buffer, data->buffer.data + data->write_ptr, bytesWritten);
		data->write_ptr += bytesWritten;
		data->buffer.length -= bytesWritten;
	}

	return bytesWritten;
}

void s3::WriteObject(const std::string& objectName, const buffer& buff)
{
	Initialize();

	const auto bucketAndFileName = GetBucketAndFileName(objectName);
	const auto bucket = bucketAndFileName[0];
	const auto key = bucketAndFileName[1];

	const char* host = getenv("S3_HOSTNAME");

	if (!host)
	{
		throw std::runtime_error("Environment variable S3_HOSTNAME not defined");
	}

#ifdef S3_DEFAULT_REGION

	// extract region name from host name, assuming aws
	// s3.us-east-1.amazonaws.com
	auto tokens = util::Split(host, ".");

	if (tokens.size() == 3)
	{
		std::cerr << "Hostname does not follow pattern s3.<regionname>.amazonaws.com" << std::endl;
		return;
	}

	const std::string region = tokens[1];

	// libs3 boilerplate

	// clang-format off

        S3BucketContext bucketContext =
        {
                host,
                bucket.c_str(),
                S3ProtocolHTTP,
                S3UriStylePath,
                access_key,
                secret_key,
                security_token,
                region.c_str()
        };
#else

	// clang-format off

	S3BucketContext bucketContext =
	{
		host,
		bucket.c_str(),
		S3ProtocolHTTP,
		S3UriStylePath,
		access_key,
		secret_key,
		security_token
	};
#endif

	// clang-format on

	S3PutObjectHandler putObjectHandler = {responseHandler, &putObjectDataCallback};

	write_data data;
	data.buffer.data = buff.data;
	data.buffer.length = buff.length;
	data.write_ptr = 0;

	logger logr("s3");

	timer t(true);

	int count = 0;
	do
	{
		if (count > 0)
		{
			sleep(2 * count);
		}

		logr.Debug("Writing to host=" + std::string(host) + " bucket=" + bucket + " key=" + key + " (" +
		           S3_get_status_name(statusG) + ")");

#ifdef S3_DEFAULT_REGION
		S3_put_object(&bucketContext, key.c_str(), buff.length, NULL, NULL, 0, &putObjectHandler, &data);
#else
		S3_put_object(&bucketContext, key.c_str(), buff.length, NULL, NULL, &putObjectHandler, &data);
#endif

		count++;
	} while (S3_status_is_retryable(statusG) && count < 3);

	// remove pointer to original buff so that double free doesn't occur
	data.buffer.data = 0;

	switch (statusG)
	{
		case S3StatusOK:
		{
			t.Stop();
			const double time = static_cast<double>(t.GetTime());
			const double size = util::round(static_cast<double>(buff.length) / 1024. / 1024., 1);
			logr.Info("Wrote " + std::to_string(size) + "MB in " + std::to_string(time) + " ms (" +
			          std::to_string(size / time) + " MB/s)");
		}
		break;
		case S3StatusInternalError:
			logr.Error(std::string(S3_get_status_name(statusG)) + ": is there a proxy blocking the connection?");
			throw himan::kFileDataNotFound;
		case S3StatusFailedToConnect:
			logr.Error(std::string(S3_get_status_name(statusG)) + ": is proxy required but not set?");
			throw himan::kFileDataNotFound;
		case S3StatusErrorInvalidAccessKeyId:
			logr.Error(std::string(S3_get_status_name(statusG)) +
			           ": are Temporary Security Credentials used without security token (env: S3_SESSION_TOKEN)?");
			throw himan::kFileDataNotFound;
		default:
			logr.Error(S3_get_status_name(statusG));
			throw himan::kFileDataNotFound;
	}
}

#else
buffer s3::ReadFile(const file_information& fileInformation)
{
	throw std::runtime_error("S3 support not compiled");
}
void s3::WriteObject(const std::string& objectName, const buffer& buff)
{
	throw std::runtime_error("S3 support not compiled");
}
#endif
