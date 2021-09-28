#include "drought.h"
#include "commonConstants.h"
#include "basicMath.h"
#include "gammaFunction.h"

Drought::Drought(droughtIndex index, int firstYear, int lastYear, Crit3DDate date, Crit3DMeteoPoint* meteoPoint, Crit3DMeteoSettings* meteoSettings)
{
    this->index = index;
    this->firstYear = firstYear;
    this->lastYear = lastYear;
    this->date = date;
    this->meteoPoint = meteoPoint;
    this->meteoSettings = meteoSettings;
    this->timeScale = 3; //default
    this->computeAll = false;  //default
}

droughtIndex Drought::getIndex() const
{
    return index;
}

void Drought::setIndex(const droughtIndex &value)
{
    index = value;
}

int Drought::getTimeScale() const
{
    return timeScale;
}

void Drought::setTimeScale(int value)
{
    timeScale = value;
}

int Drought::getFirstYear() const
{
    return firstYear;
}

void Drought::setFirstYear(int value)
{
    firstYear = value;
}

int Drought::getLastYear() const
{
    return lastYear;
}

void Drought::setLastYear(int value)
{
    lastYear = value;
}

bool Drought::getComputeAll() const
{
    return computeAll;
}

void Drought::setComputeAll(bool value)
{
    computeAll = value;
}

float Drought::computeDroughtIndex()
{
    timeScale = timeScale - 1; // index start from 0
    if (index == INDEX_SPI)
    {
        if (!computeSpiParameters())
        {
            return NODATA;
        }
    }
    else if (index == INDEX_SPEI)
    {
        if (!computeSpeiParameters())
        {
            return NODATA;
        }
    }
    int start;
    int end;
    std::vector<float> mySum(meteoPoint->nrObsDataDaysM);
    std::vector<float> myResults(meteoPoint->nrObsDataDaysM);
    if (computeAll)
    {
        start = timeScale;
        end = meteoPoint->nrObsDataDaysM;
        for (int i = 0; i <= timeScale; i++)
        {
            mySum[i] = NODATA;
            myResults[i] = NODATA;
        }
    }
    else
    {
        int currentYear = date.year;
        int currentMonth = date.month;
        end = (currentYear - meteoPoint->obsDataM[0]._year)*12 + currentMonth-meteoPoint->obsDataM[0]._month + 1; // come nel caso precedente end al massimo è pari a meteoPoint->nrObsDataDaysM
        start = end - 1;
        if (end > meteoPoint->nrObsDataDaysM)
        {
            return NODATA;
        }
        for (int i = 0; i <= start-timeScale; i++)
        {
            mySum[i] = NODATA;
            myResults[i] = NODATA;
        }
    }

    for (int j = start; j < end; j++)
    {
        mySum[j] = 0;
        for(int i = 0; i<=timeScale; i++)
        {
            if ((j-i)>=0 && j < meteoPoint->nrObsDataDaysM)
            {
                if (index == INDEX_SPI)
                {
                    if (meteoPoint->obsDataM[j-i].prec != NODATA)
                    {
                        mySum[j] = mySum[j] + meteoPoint->obsDataM[j-i].prec;
                    }
                    else
                    {
                        mySum[j] = NODATA;
                        break;
                    }
                }
                else if(index == INDEX_SPEI)
                {
                    if (meteoPoint->obsDataM[j-i].prec != NODATA && meteoPoint->obsDataM[j-i].et0_hs != NODATA)
                    {
                        mySum[j] = mySum[j] + meteoPoint->obsDataM[j-i].prec - meteoPoint->obsDataM[j-i].et0_hs;
                    }
                    else
                    {
                        mySum[j] = NODATA;
                        break;
                    }
                }
            }
            else
            {
                mySum[j] = NODATA;
            }
        }
    }

    for (int j = start; j < end; j++)
    {
        int myMonthIndex = (j % 12)+1;  //start from 1
        myResults[j] = NODATA;

        if (mySum[j] != NODATA)
        {
            if (index == INDEX_SPI)
            {
                // TO DO
                // myResults[j] = standardGaussianInvCDF(math.gammaCDF(mySum(j), currentGamma(myMonthIndex))) // già presente in math
            }
            else if(index == INDEX_SPEI)
            {
                // TO DO
                // myResults[j] = standardGaussianInvCDF(math.logLogisticCDF(mySum(j), currentLogLogistic(myMonthIndex))) // già presente in math
            }
        }
    }
    return myResults[end-1];
}

bool Drought::computeSpiParameters()
{
    if (meteoPoint->nrObsDataDaysM == 0)
    {
        return false;
    }

    if (meteoPoint->obsDataM[0]._year > lastYear || meteoPoint->obsDataM[meteoPoint->nrObsDataDaysM-1]._year < firstYear)
    {
        return false;
    }
    int indexStart;
    if (firstYear == meteoPoint->obsDataM[0]._year)
    {
        indexStart = timeScale;
    }
    else
    {
        indexStart = (firstYear - meteoPoint->obsDataM[0]._year)*12 - (meteoPoint->obsDataM[0]._month-1);
        if (indexStart < timeScale)
        {
            indexStart = timeScale;
        }
    }
    if (meteoPoint->obsDataM[indexStart]._year > lastYear)
    {
        return false;
    }

    // int firstYearStation = std::max(meteoPoint->obsDataM[indexStart]._year, firstYear); // LC non viene mai usata nel codice vb
    int lastYearStation = std::min(meteoPoint->obsDataM[meteoPoint->nrObsDataDaysM-1]._year, lastYear);

    int n = 0;
    float count = 0;
    int nTot = 0;
    std::vector<float> mySums;
    std::vector<float> monthSeries;

    for (int j = indexStart; j<meteoPoint->nrObsDataDaysM; j++)
    {
        if (meteoPoint->obsDataM[j]._year <= lastYearStation)
        {
            count = 0;
            nTot = 0;
            mySums.push_back(0);
            for(int i = 0; i<= timeScale; i++)
            {
                nTot = nTot + 1;
                if (meteoPoint->obsDataM[j-i].prec != NODATA)
                {
                    mySums[n] = mySums[n] + meteoPoint->obsDataM[j-i].prec;
                    count = count + 1;
                }
                else
                {
                        mySums[n] = NODATA;
                        count = 0;
                        break;
                }
            }
            if ( (count / nTot) < (meteoSettings->getMinimumPercentage() / 100) )
            {
                mySums[n] = NODATA;
            }
            n = n + 1;
        }
        else
        {
            break;
        }
    }

    for (int i = 0; i<12; i++)
    {
        int myMonth = (meteoPoint->obsDataM[indexStart]._month + i) % 12;  //start from 1

        gammaStruct.beta = NODATA;
        gammaStruct.gamma = NODATA;
        gammaStruct.pzero = NODATA;
        currentGamma.push_back(gammaStruct);

        for (int j=i; j<mySums.size(); j=j+12)
        {
            if (mySums[j] != NODATA)
            {
                monthSeries.push_back(mySums[j]);
            }
        }

        if (monthSeries.size() / (mySums.size()/12) >= meteoSettings->getMinimumPercentage() / 100)
        {
            gammaDistributions::gammaFitting(monthSeries, n, &(currentGamma[myMonth].beta), &(currentGamma[myMonth].gamma),  &(currentGamma[myMonth].pzero));
        }
    }
    return true;
}

bool Drought::computeSpeiParameters()
{
    if (meteoPoint->nrObsDataDaysM == 0)
    {
        return false;
    }

    if (meteoPoint->obsDataM[0]._year > lastYear || meteoPoint->obsDataM[meteoPoint->nrObsDataDaysM-1]._year < firstYear)
    {
        return false;
    }
    int indexStart;
    if (firstYear == meteoPoint->obsDataM[0]._year)
    {
        indexStart = timeScale;
    }
    else
    {
        indexStart = (firstYear - meteoPoint->obsDataM[0]._year)*12 - (meteoPoint->obsDataM[0]._month-1);
        if (indexStart < timeScale)
        {
            indexStart = timeScale;
        }
    }
    if (meteoPoint->obsDataM[indexStart]._year > lastYear)
    {
        return false;
    }

    // int firstYearStation = std::max(meteoPoint->obsDataM[indexStart]._year, firstYear); // LC non viene mai usata nel codice vb
    int lastYearStation = std::min(meteoPoint->obsDataM[meteoPoint->nrObsDataDaysM-1]._year, lastYear);

    int n = 0;
    float count = 0;
    int nTot = 0;
    std::vector<float> mySums;
    std::vector<float> monthSeries;

    for (int j = indexStart; j<meteoPoint->nrObsDataDaysM; j++)
    {
        if (meteoPoint->obsDataM[j]._year <= lastYearStation)
        {
            count = 0;
            nTot = 0;
            mySums.push_back(0);
            for(int i = 0; i<=timeScale; i++)
            {
                nTot = nTot + 1;
                if (meteoPoint->obsDataM[j-i].prec != NODATA && meteoPoint->obsDataM[j-i].et0_hs != NODATA)
                {
                    mySums[n] = mySums[n] + meteoPoint->obsDataM[j-i].prec - meteoPoint->obsDataM[j-i].et0_hs;
                    count = count + 1;
                }
                else
                {
                        mySums[n] = NODATA;
                        count = 0;
                        break;
                }
            }
            if ( (count / nTot) < (meteoSettings->getMinimumPercentage() / 100))
            {
                mySums[n] = NODATA;
            }
            n = n + 1;
        }
        else
        {
            break;
        }
    }

    for (int i = 0; i<12; i++)
    {
        int myMonth = (meteoPoint->obsDataM[indexStart]._month + i) % 12; //start from 1

        logLogisticStruct.beta = NODATA;
        logLogisticStruct.gamma = NODATA;
        logLogisticStruct.alpha = NODATA;
        currentLogLogistic.push_back(logLogisticStruct);

        for (int j=i; j<mySums.size(); j=j+12)
        {
            if (mySums[j] != NODATA)
            {
                monthSeries.push_back(mySums[j]);
            }
        }

        if (monthSeries.size() / (mySums.size()/12) >= meteoSettings->getMinimumPercentage() / 100)
        {
            // TO DO
                // Sort values
                sorting::quicksortAscendingFloat(monthSeries, 0, unsigned(monthSeries.size() - 1));
                /*
                // Compute probability weighted moments
                math.probabilityWeightedMoments monthSeries, n, PWM, 0, 0, False
                // Fit a Log Logistic probability function
                math.logLogisticFitting PWM, currentLogLogistic(myMonth)
            */
        }
    }
    return true;
}

void Drought::setMeteoPoint(Crit3DMeteoPoint *value)
{
    meteoPoint = value;
}

Crit3DMeteoSettings *Drought::getMeteoSettings() const
{
    return meteoSettings;
}

Crit3DDate Drought::getDate() const
{
    return date;
}

void Drought::setDate(const Crit3DDate &value)
{
    date = value;
}
