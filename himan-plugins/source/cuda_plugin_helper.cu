#include "cuda_plugin_helper.h"
#include "plugin_factory.h"
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/system/cuda/execution_policy.h>

#define HIMAN_AUXILIARY_INCLUDE

#include "cache.h"
#include "fetcher.h"

#undef HIMAN_AUXILIARY_INCLUDE

namespace himan
{
namespace cuda
{
bool Unpack(std::shared_ptr<himan::info<double>> fullInfo, cudaStream_t& stream, double* d_arr)
{
	return Unpack<double>(fullInfo, stream, d_arr);
}

template <typename T>
bool Unpack(std::shared_ptr<himan::info<T>> fullInfo, cudaStream_t& stream, T* d_arr)
{
	using namespace himan;
	using namespace himan::plugin;

	const size_t N = fullInfo->SizeLocations();

	ASSERT(N > 0);
	ASSERT(d_arr);

	// Unpack if needed, leave data to device
	if (fullInfo->PackedData()->HasData())
	{
		ASSERT(fullInfo->PackedData()->packingType == kSimplePacking);
		ASSERT(fullInfo->Data().Size() == N);

		packing::Unpack(std::dynamic_pointer_cast<simple_packed>(fullInfo->PackedData()).get(), d_arr, &stream);

		return true;
	}
	else
	{
		// Data was not packed, ie it was returned to us from cache
		CUDA_CHECK(
		    cudaMemcpyAsync(d_arr, fullInfo->Data().ValuesAsPOD(), sizeof(T) * N, cudaMemcpyHostToDevice, stream));

		CUDA_CHECK(cudaStreamSynchronize(stream));
		return false;
	}
}

template bool Unpack<double>(std::shared_ptr<himan::info<double>>, cudaStream_t&, double*);
template bool Unpack<float>(std::shared_ptr<himan::info<float>>, cudaStream_t&, float*);

template <typename T>
void PrepareInfo(std::shared_ptr<himan::info<T>> info, T* d_ret, cudaStream_t& stream, bool copyToHost)
{
	if (Unpack(info, stream, d_ret) && copyToHost)
	{
		CUDA_CHECK(cudaMemcpyAsync(info->Data().ValuesAsPOD(), d_ret, sizeof(T) * info->SizeLocations(),
		                           cudaMemcpyDeviceToHost, stream));
		CUDA_CHECK(cudaStreamSynchronize(stream));

		info->PackedData()->Clear();

		auto c = GET_PLUGIN(cache);
		c->Insert(info);
	}
}

template void PrepareInfo<double>(std::shared_ptr<himan::info<double>>, double*, cudaStream_t&, bool);
template void PrepareInfo<float>(std::shared_ptr<himan::info<float>>, float*, cudaStream_t&, bool);

/*
template <>
void PrepareInfo(std::shared_ptr<himan::info<double>> info, float* d_ret, cudaStream_t& stream, bool copyToHost)
{
    const size_t N = info->SizeLocations();
    double* d_arr = 0;
    CUDA_CHECK(cudaMalloc(reinterpret_cast<double**>(&d_arr), N * sizeof(double)));

    if (Unpack(info, stream, d_arr) && copyToHost)
    {
        CUDA_CHECK(cudaMemcpyAsync(info->Data().ValuesAsPOD(), d_arr, sizeof(double) * info->SizeLocations(),
                                   cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));

        info->PackedData()->Clear();

        auto c = GET_PLUGIN(cache);
        c->Insert(info);
    }

    thrust::device_ptr<double> dt_arr = thrust::device_pointer_cast(d_arr);
    thrust::device_ptr<float> dt_farr = thrust::device_pointer_cast(d_ret);

    thrust::copy(thrust::cuda::par.on(stream), dt_arr, dt_arr + N, dt_farr);
    thrust::replace_if(thrust::cuda::par.on(stream), dt_farr, dt_farr + N,
                       [] __device__(const float& val) { return ::isnan(val); }, himan::MissingFloat());

    CUDA_CHECK(cudaStreamSynchronize(stream));
    CUDA_CHECK(cudaFree(d_arr));
}
*/

template <typename T>
void ReleaseInfo(std::shared_ptr<himan::info<T>> info, T* d_arr, cudaStream_t& stream)
{
	CUDA_CHECK(cudaMemcpyAsync(info->Data().ValuesAsPOD(), d_arr, info->SizeLocations() * sizeof(T),
	                           cudaMemcpyDeviceToHost, stream));
	CUDA_CHECK(cudaStreamSynchronize(stream));
}

template void ReleaseInfo<double>(std::shared_ptr<himan::info<double>>, double*, cudaStream_t&);
template void ReleaseInfo<float>(std::shared_ptr<himan::info<float>>, float*, cudaStream_t&);

/*
template <>
void ReleaseInfo(std::shared_ptr<himan::info<double>> info, float* d_arr, cudaStream_t& stream)
{
    const size_t N = info->SizeLocations();

    float* h_arr = new float[N];
    CUDA_CHECK(cudaMemcpyAsync(h_arr, d_arr, info->SizeLocations() * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));

    auto& res = VEC(info);

    std::copy(h_arr, h_arr + N, res.begin());
    std::replace_if(res.begin(), res.end(), [](const double& val) { return ::isnan(val); }, himan::MissingDouble());

    delete[] h_arr;
}
*/
template <typename T>
std::shared_ptr<himan::info<T>> Fetch(const std::shared_ptr<const plugin_configuration> conf,
                                      const himan::forecast_time& theTime, const himan::level& theLevel,
                                      const himan::params& theParams, const himan::forecast_type& theType,
                                      bool returnPacked)
{
	for (const auto& p : theParams)
	{
		auto ret = Fetch<T>(conf, theTime, theLevel, p, theType, returnPacked);

		if (ret)
		{
			return ret;
		}
	}
	return nullptr;
}

template std::shared_ptr<himan::info<double>> Fetch<double>(const std::shared_ptr<const plugin_configuration>,
                                                            const himan::forecast_time&, const himan::level&,
                                                            const himan::params&, const himan::forecast_type&, bool);
template std::shared_ptr<himan::info<float>> Fetch<float>(const std::shared_ptr<const plugin_configuration>,
                                                          const himan::forecast_time&, const himan::level&,
                                                          const himan::params&, const himan::forecast_type&, bool);

template <typename T>
std::shared_ptr<himan::info<T>> Fetch(const std::shared_ptr<const plugin_configuration> conf,
                                      const himan::forecast_time& theTime, const himan::level& theLevel,
                                      const himan::param& theParam, const himan::forecast_type& theType,
                                      bool returnPacked)
{
	try
	{
		auto f = GET_PLUGIN(fetcher);
		return f->Fetch<T>(conf, theTime, theLevel, theParam, theType, returnPacked);
	}
	catch (HPExceptionType& e)
	{
		if (e != kFileDataNotFound)
		{
			throw std::runtime_error("cape_cuda::Fetch(): Unable to proceed");
		}

		return nullptr;
	}
}

template std::shared_ptr<himan::info<double>> Fetch<double>(const std::shared_ptr<const plugin_configuration>,
                                                            const himan::forecast_time&, const himan::level&,
                                                            const himan::param&, const himan::forecast_type&, bool);
template std::shared_ptr<himan::info<float>> Fetch<float>(const std::shared_ptr<const plugin_configuration>,
                                                          const himan::forecast_time&, const himan::level&,
                                                          const himan::param&, const himan::forecast_type&, bool);

}  // namespace cuda
}  // namespace himan
