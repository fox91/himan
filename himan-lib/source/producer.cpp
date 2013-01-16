/*
 * producer.cpp
 *
 *  Created on: Dec 28, 2012
 *      Author: partio
 */

#include "producer.h"

using namespace himan;

producer::producer()
    : itsFmiProducerId(kHPMissingInt)
    , itsProcess(kHPMissingInt)
    , itsCentre(kHPMissingInt)
    , itsNeonsName("")
{

}

producer::producer(long theFmiProducerId)
    : itsFmiProducerId(theFmiProducerId)
    , itsProcess(kHPMissingInt)
    , itsCentre(kHPMissingInt)
    , itsNeonsName("")
{

}

producer::producer(long theCentre, long theProcess)
    : itsFmiProducerId(kHPMissingInt)
    , itsProcess(theProcess)
    , itsCentre(theCentre)
    , itsNeonsName("")
{

}

producer::producer(long theFmiProducerId, long theCentre, long theProcess, const std::string& theNeonsName)
    : itsFmiProducerId(theFmiProducerId)
    , itsProcess(theProcess)
    , itsCentre(theCentre)
    , itsNeonsName("")
{

}

void producer::Centre(long theCentre)
{
    itsCentre = theCentre;
}

long producer::Centre() const
{
    return itsCentre;
}

void producer::Process(long theProcess)
{
    itsProcess = theProcess;
}

long producer::Process() const
{
    return itsProcess;
}

void producer::Id(long theId)
{
    itsFmiProducerId = theId;
}

long producer::Id() const
{
    return itsFmiProducerId;
}

void producer::Name(const std::string& theName)
{
    itsNeonsName = theName;
}

std::string producer::Name() const
{
    return itsNeonsName;
}

std::ostream& producer::Write(std::ostream& file) const
{

    file << "<" << ClassName() << " " << Version() << ">" << std::endl;

    file << "__itsFmiProducerId__ " << itsFmiProducerId << std::endl;
    file << "__itsProcess__ " << itsProcess << std::endl;
    file << "__itsCentre__ " << itsCentre << std::endl;
    file << "__itsNeonsName__ " << itsNeonsName << std::endl;

    return file;

}
