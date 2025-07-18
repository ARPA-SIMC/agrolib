#include <math.h>

#include "commonConstants.h"
#include "basicMath.h"
#include "spatialControl.h"
#include "interpolation.h"
#include "statistics.h"

float findThreshold(meteoVariable myVar, Crit3DMeteoSettings* meteoSettings,
                    float value, float stdDev, float nrStdDev, float avgDeltaZ, float minDistance)
{
    float zWeight, distWeight, threshold;

    if (   myVar == precipitation
        || myVar == dailyPrecipitation)
    {
        distWeight = MAXVALUE(1.f, minDistance / 2000.f);
        if (value <= meteoSettings->getRainfallThreshold())
            threshold = MAXVALUE(5.f, distWeight + stdDev * (nrStdDev + 1));
        else
            return 900.f;
    }
    else if (   myVar == airTemperature
             || myVar == airDewTemperature
             || myVar == dailyAirTemperatureMax
             || myVar == dailyAirTemperatureMin
             || myVar == dailyAirTemperatureAvg )
    {
        threshold = 1.f;
        zWeight = avgDeltaZ / 100.f;
        distWeight = minDistance / 1000.f;

        threshold = MINVALUE(MINVALUE(distWeight + threshold + zWeight, 12.f) + stdDev * nrStdDev, 15.f);
    }
    else if (   myVar == airRelHumidity
             || myVar == dailyAirRelHumidityMax
             || myVar == dailyAirRelHumidityMin
             || myVar == dailyAirRelHumidityAvg )
    {
        threshold = 20.f;
        zWeight = avgDeltaZ / 10.f;
        distWeight = minDistance / 1000.f;
        threshold += zWeight + distWeight + stdDev * nrStdDev;
    }
    else if (   myVar == windScalarIntensity
             || myVar == windVectorIntensity
             || myVar == dailyWindScalarIntensityAvg
             || myVar == dailyWindScalarIntensityMax
             || myVar == dailyWindVectorIntensityAvg
             || myVar == dailyWindVectorIntensityMax)
    {
        threshold = 1.f;
        zWeight = avgDeltaZ / 50.f;
        distWeight = minDistance / 2000.f;
        threshold += zWeight + distWeight + stdDev * nrStdDev;
    }
    else if (   myVar == globalIrradiance)
    {
        threshold = 500;
        distWeight = minDistance / 5000.f;
        threshold += distWeight + stdDev * (nrStdDev + 1.f);
    }
    else if (   myVar ==  dailyGlobalRadiation)
    {
        threshold = 10;
        distWeight = minDistance / 5000.f;
        threshold += distWeight + stdDev * (nrStdDev + 1.f);
    }
    else if (myVar == atmTransmissivity)
        threshold = MAXVALUE(stdDev * nrStdDev, 0.25f);
    else
        threshold = stdDev * nrStdDev;

    return threshold;
}


bool computeResiduals(meteoVariable myVar, Crit3DMeteoPoint* meteoPoints, int nrMeteoPoints,
                      std::vector <Crit3DInterpolationDataPoint> &interpolationPoints, Crit3DInterpolationSettings* settings,
                      Crit3DMeteoSettings* meteoSettings, bool excludeOutsideDem, bool excludeSupplemental)
{

    if (myVar == noMeteoVar) return false;

    std::vector <double> myProxyValues;
    bool isValid;

    for (int i = 0; i < nrMeteoPoints; i++)
    {
        myProxyValues = meteoPoints[i].getProxyValues();

        meteoPoints[i].residual = NODATA;

        if (meteoPoints[i].active)
        {
            isValid = (! excludeSupplemental || checkLapseRateCode(meteoPoints[i].lapseRateCode, settings->getUseLapseRateCode(), false));
            isValid = (isValid && (! excludeOutsideDem || meteoPoints[i].isInsideDem));

            if (isValid && meteoPoints[i].quality == quality::accepted)
            {
                float myValue = meteoPoints[i].currentValue;

                float interpolatedValue = interpolate(interpolationPoints, settings, meteoSettings, myVar,
                                                      float(meteoPoints[i].point.utm.x),
                                                      float(meteoPoints[i].point.utm.y),
                                                      float(meteoPoints[i].point.z),
                                                      myProxyValues, false);

                if (  myVar == precipitation || myVar == dailyPrecipitation)
                {
                    if (myValue != NODATA)
                    {
                        if (myValue < meteoSettings->getRainfallThreshold())
                            myValue=0.;
                    }

                    if (interpolatedValue != NODATA)
                    {
                        if (interpolatedValue < meteoSettings->getRainfallThreshold())
                            interpolatedValue=0.;
                    }
                }

                // TODO derived var

                if ((interpolatedValue != NODATA) && (myValue != NODATA))
                {
                    meteoPoints[i].residual = myValue - interpolatedValue;
                }
            }
        }
    }

    return true;
}


bool computeResidualsLocalDetrending(meteoVariable myVar, const Crit3DTime &myTime, Crit3DMeteoPoint* meteoPoints, int nrMeteoPoints,
                                     std::vector <Crit3DInterpolationDataPoint> &interpolationPoints,
                                     Crit3DInterpolationSettings* settings,
                                     Crit3DMeteoSettings* meteoSettings, Crit3DClimateParameters* climateParameters,
                                     bool excludeOutsideDem, bool excludeSupplemental)
{

    if (myVar == noMeteoVar) return false;

    std::vector <double> myProxyValues;
    bool isValid;
    std::string errorStdString;

    for (int i = 0; i < nrMeteoPoints; i++)
    {
        myProxyValues = meteoPoints[i].getProxyValues();

        meteoPoints[i].residual = NODATA;

        if (meteoPoints[i].active)
        {

            isValid = (! excludeSupplemental || checkLapseRateCode(meteoPoints[i].lapseRateCode, settings->getUseLapseRateCode(), false));
            isValid = (isValid && (! excludeOutsideDem || meteoPoints[i].isInsideDem));

            if (isValid && meteoPoints[i].quality == quality::accepted)
            {
                float myValue = meteoPoints[i].currentValue;

                std::vector <Crit3DInterpolationDataPoint> subsetInterpolationPoints;
                localSelection(interpolationPoints, subsetInterpolationPoints, float(meteoPoints[i].point.utm.x),
                               float(meteoPoints[i].point.utm.y), *settings, false);
                if (! preInterpolation(subsetInterpolationPoints, settings, meteoSettings,
                                      climateParameters, meteoPoints, nrMeteoPoints, myVar, myTime, errorStdString))
                {
                    return false;
                }

                float interpolatedValue = interpolate(subsetInterpolationPoints, settings, meteoSettings, myVar,
                                                      float(meteoPoints[i].point.utm.x),
                                                      float(meteoPoints[i].point.utm.y),
                                                      float(meteoPoints[i].point.z),
                                                      myProxyValues, false);

                if (  myVar == precipitation || myVar == dailyPrecipitation)
                {
                    if (myValue != NODATA)
                    {
                        if (myValue < meteoSettings->getRainfallThreshold())
                            myValue=0.;
                    }

                    if (interpolatedValue != NODATA)
                    {
                        if (interpolatedValue < meteoSettings->getRainfallThreshold())
                            interpolatedValue=0.;
                    }
                }

                // TODO derived var

                if ((interpolatedValue != NODATA) && (myValue != NODATA))
                {
                    meteoPoints[i].residual = myValue - interpolatedValue;
                }
            }
        }
    }

    return true;
}


bool computeResidualsGlocalDetrending(meteoVariable myVar, const Crit3DMacroArea &myArea, int elevationPos,
                                      Crit3DMeteoPoint* meteoPoints, std::vector <Crit3DInterpolationDataPoint> &interpolationPoints,
                                      Crit3DInterpolationSettings* settings, Crit3DMeteoSettings* meteoSettings,
                                      bool excludeOutsideDem, bool excludeSupplemental)
{
    //TODO: glocal cv with grid ONLY (no DEM)
    std::vector<Crit3DInterpolationDataPoint> areaInterpolationPoints;
    std::vector<int> meteoPointsList = myArea.getMeteoPoints();
    std::vector<float> areaCells;
    std::vector<double> myProxyValues;
    bool isValid;

    //un solo detrending per ogni area
    macroAreaDetrending(myArea, myVar, *settings, meteoSettings, meteoPoints, interpolationPoints, areaInterpolationPoints, elevationPos);

    //ciclo sui meteopoint dell'area
    for (int i = 0; i < (int)meteoPointsList.size(); i++)
    {
        myProxyValues = meteoPoints[meteoPointsList[i]].getProxyValues();

        if (meteoPoints[meteoPointsList[i]].active)
        {
            //peso della stazione nell'area attuale in base alla sua posizione
            int row, col;
            float weight = NODATA;
            int temp = NODATA;

            //valido solo per DEM
            areaCells = myArea.getAreaCellsDEM();
            std::string name = meteoPoints[meteoPointsList[i]].name;
            std::string id = meteoPoints[meteoPointsList[i]].id;

            gis::getRowColFromXY(*settings->getCurrentDEM()->header, meteoPoints[meteoPointsList[i]].point.utm, &row, &col);
            temp = settings->getCurrentDEM()->header->nrCols*row + col;

            for (int k = 0; k < (int)areaCells.size(); k = k + 2)
            {
                if (areaCells[k] == temp)
                    weight = areaCells[k+1];
            }

            isValid = (! excludeSupplemental || checkLapseRateCode(meteoPoints[meteoPointsList[i]].lapseRateCode, settings->getUseLapseRateCode(), false));
            isValid = (isValid && (! excludeOutsideDem || meteoPoints[meteoPointsList[i]].isInsideDem));

            if (isValid && meteoPoints[meteoPointsList[i]].quality == quality::accepted)
            {
                float myValue = meteoPoints[meteoPointsList[i]].currentValue;

                float interpolatedValue = interpolate(areaInterpolationPoints, settings, meteoSettings, myVar,
                                                      float(meteoPoints[meteoPointsList[i]].point.utm.x),
                                                      float(meteoPoints[meteoPointsList[i]].point.utm.y),
                                                      float(meteoPoints[meteoPointsList[i]].point.z),
                                                      myProxyValues, false);


                if (!isEqual(interpolatedValue, NODATA) && !isEqual(myValue, NODATA) && !isEqual(weight, NODATA))
                {
                    if (isEqual(meteoPoints[meteoPointsList[i]].residual, NODATA))
                        meteoPoints[meteoPointsList[i]].residual = (myValue - interpolatedValue)*weight;
                    else
                    {
                        meteoPoints[meteoPointsList[i]].residual += (myValue - interpolatedValue)*weight;
                    }

                }
            }
        }
    }

    return true;
}


float computeErrorCrossValidation(Crit3DMeteoPoint* myPoints, int nrMeteoPoints)
{
    std::vector <float> obsValues, estValues;

    for (int i=0; i < nrMeteoPoints; i++)
    {
        if (myPoints[i].active)
        {
            float value = myPoints[i].currentValue;
            float residual = myPoints[i].residual;

            if (value != NODATA && residual != NODATA)
            {
                obsValues.push_back(value);
                estValues.push_back(value - residual);
            }
        }
    }

    if (obsValues.size() > 0)
    {
        return statistics::meanAbsoluteError(obsValues, estValues);
    }
    else
        return NODATA;
}


bool spatialQualityControl(meteoVariable myVar, Crit3DMeteoPoint* meteoPoints, int nrMeteoPoints,
                           Crit3DInterpolationSettings *settings, Crit3DMeteoSettings* meteoSettings,
                           Crit3DClimateParameters* myClimate, const Crit3DTime &myTime, std::string &errorStr)
{
    float stdDev, avgDeltaZ, minDist, myValue, myResidual;
    std::vector <int> listIndex;
    std::vector <float> listResiduals;
    std::vector <Crit3DInterpolationDataPoint> myInterpolationPoints;

    if (passDataToInterpolation(meteoPoints, nrMeteoPoints, myInterpolationPoints, settings))
    {
        // detrend
        if (! preInterpolation(myInterpolationPoints, settings, meteoSettings, myClimate,
                              meteoPoints, nrMeteoPoints, myVar, myTime, errorStr))
        {
            return false;
        }

        // compute residuals
        if (! computeResiduals(myVar, meteoPoints, nrMeteoPoints, myInterpolationPoints, settings, meteoSettings, false, false))
        {
            errorStr = "Error in compute residuals.";
            return false;
        }

        int i;
        for (i = 0; i < nrMeteoPoints; i++)
        {
            if (meteoPoints[i].quality == quality::accepted)
            {
                if (neighbourhoodVariability(myVar, myInterpolationPoints, settings, float(meteoPoints[i].point.utm.x),
                         float(meteoPoints[i].point.utm.y),float(meteoPoints[i].point.z),
                         10, &stdDev, &avgDeltaZ, &minDist))
                {
                    myValue = meteoPoints[i].currentValue;
                    myResidual = meteoPoints[i].residual;
                    stdDev = MAXVALUE(stdDev, myValue/100.f);
                    if (fabs(myResidual) > findThreshold(myVar, meteoSettings, myValue, stdDev, 2, avgDeltaZ, minDist))
                    {
                        listIndex.push_back(i);
                        meteoPoints[i].quality = quality::wrong_spatial;
                    }
                }
            }
        }

        if (listIndex.size() > 0)
        {
            if (passDataToInterpolation(meteoPoints, nrMeteoPoints, myInterpolationPoints, settings))
            {
                if (! preInterpolation(myInterpolationPoints, settings, meteoSettings, myClimate,
                                      meteoPoints, nrMeteoPoints, myVar, myTime, errorStr))
                {
                    return false;
                }

                float interpolatedValue;
                for (i=0; i < int(listIndex.size()); i++)
                {
                    interpolatedValue = interpolate(myInterpolationPoints, settings, meteoSettings, myVar,
                                            float(meteoPoints[listIndex[i]].point.utm.x),
                                            float(meteoPoints[listIndex[i]].point.utm.y),
                                            float(meteoPoints[listIndex[i]].point.z),
                                            meteoPoints[i].getProxyValues(),
                                            false);

                    myValue = meteoPoints[listIndex[i]].currentValue;

                    listResiduals.push_back(interpolatedValue - myValue);
                }

                for (i=0; i < int(listIndex.size()); i++)
                {
                    if (neighbourhoodVariability(myVar, myInterpolationPoints, settings, float(meteoPoints[listIndex[i]].point.utm.x),
                             float(meteoPoints[listIndex[i]].point.utm.y),
                             float(meteoPoints[listIndex[i]].point.z),
                             10, &stdDev, &avgDeltaZ, &minDist))
                    {
                        myResidual = listResiduals[i];

                        myValue = meteoPoints[listIndex[i]].currentValue;

                        if (fabs(myResidual) > findThreshold(myVar, meteoSettings, myValue, stdDev, 3, avgDeltaZ, minDist))
                            meteoPoints[listIndex[i]].quality = quality::wrong_spatial;
                        else
                            meteoPoints[listIndex[i]].quality = quality::accepted;
                    }
                    else
                        meteoPoints[listIndex[i]].quality = quality::accepted;
                }
            }
        }
    }

    return true;
}


bool checkData(Crit3DQuality* myQuality, meteoVariable myVar, Crit3DMeteoPoint* meteoPoints, int nrMeteoPoints,
               const Crit3DTime &myTime, Crit3DInterpolationSettings* spatialQualityInterpolationSettings,
               Crit3DMeteoSettings* meteoSettings, Crit3DClimateParameters* myClimate, bool checkSpatial, std::string &errorStr)
{
    if (nrMeteoPoints == 0)
        return false;

    if (myVar == elaboration)
    {
        // assign data
        for (int i = 0; i < nrMeteoPoints; i++)
        {
            meteoPoints[i].currentValue = meteoPoints[i].elaboration;
            if (! isEqual(meteoPoints[i].currentValue, NODATA))
                meteoPoints[i].quality = quality::accepted;
            else
                meteoPoints[i].quality = quality::missing_data;
        }
    }
    else if (myVar == anomaly)
    {
        // assign data
        for (int i = 0; i < nrMeteoPoints; i++)
        {
            meteoPoints[i].currentValue = meteoPoints[i].anomaly;
            if (! isEqual(meteoPoints[i].currentValue, NODATA))
                meteoPoints[i].quality = quality::accepted;
            else
                meteoPoints[i].quality = quality::missing_data;
        }
    }
    else
    {
        // assign data
        for (int i = 0; i < nrMeteoPoints; i++)
            meteoPoints[i].currentValue = meteoPoints[i].getMeteoPointValue(myTime, myVar, meteoSettings);

        // quality control - syntactic
        myQuality->syntacticQualityControl(myVar, meteoPoints, nrMeteoPoints);

        // quality control - spatial
        if (checkSpatial && myVar != precipitation && myVar != dailyPrecipitation
                         && myVar != windVectorX && myVar != windVectorY
                         && myVar != windVectorDirection && myVar != dailyWindVectorDirectionPrevailing)
        {
            if (! spatialQualityControl(myVar, meteoPoints, nrMeteoPoints, spatialQualityInterpolationSettings,
                                       meteoSettings, myClimate, myTime, errorStr))
            {
                return false;
            }
        }
    }

    return true;
}


// check quality and pass good data to interpolation
bool checkAndPassDataToInterpolation(Crit3DQuality* myQuality, meteoVariable myVar, Crit3DMeteoPoint* meteoPoints,
                                     int nrMeteoPoints, const Crit3DTime &myTime, Crit3DInterpolationSettings *SQinterpolationSettings,
                                     Crit3DInterpolationSettings* interpolationSettings, Crit3DMeteoSettings *meteoSettings,
                                     Crit3DClimateParameters *myClimate, std::vector<Crit3DInterpolationDataPoint> &interpolationPoints,
                                     bool checkSpatial, std::string &errorStr)
{
    if (! checkData(myQuality, myVar, meteoPoints, nrMeteoPoints, myTime, SQinterpolationSettings,
                   meteoSettings, myClimate, checkSpatial, errorStr))
    {
        return false;
    }

    // return true if at least one valid data
    return passDataToInterpolation(meteoPoints, nrMeteoPoints, interpolationPoints, interpolationSettings);
}


bool passDataToInterpolation(Crit3DMeteoPoint* meteoPoints, int nrMeteoPoints,
                            std::vector <Crit3DInterpolationDataPoint> &interpolationPoints,
                            Crit3DInterpolationSettings* mySettings)
{
    float xMin, xMax, yMin, yMax;
    float valueMin, valueMax;

    bool isSelection = isSelectionPointsActive(meteoPoints, nrMeteoPoints);
    bool isFirst = true;
    int nrValid = 0;

    interpolationPoints.clear();

    for (int i = 0; i < nrMeteoPoints; i++)
    {
        if (meteoPoints[i].active && meteoPoints[i].quality == quality::accepted && (! isSelection || meteoPoints[i].selected))
        {
            Crit3DInterpolationDataPoint myPoint;

            myPoint.index = i;
            myPoint.value = meteoPoints[i].currentValue;
            myPoint.point->utm.x = meteoPoints[i].point.utm.x;
            myPoint.point->utm.y = meteoPoints[i].point.utm.y;
            myPoint.point->z = meteoPoints[i].point.z;
            myPoint.lapseRateCode = meteoPoints[i].lapseRateCode;
            myPoint.proxyValues = meteoPoints[i].proxyValues;
            myPoint.topographicDistance = meteoPoints[i].topographicDistance;
            myPoint.isActive = true;
            myPoint.isMarked = meteoPoints[i].marked;

            if (isFirst)
            {
                xMin = float(myPoint.point->utm.x);
                xMax = float(myPoint.point->utm.x);
                yMin = float(myPoint.point->utm.y);
                yMax = float(myPoint.point->utm.y);
                valueMin = myPoint.value;
                valueMax = myPoint.value;
                isFirst = false;
            }
            else
            {
                xMin = MINVALUE(xMin, (float)myPoint.point->utm.x);
                xMax = MAXVALUE(xMax, (float)myPoint.point->utm.x);
                yMin = MINVALUE(yMin, (float)myPoint.point->utm.y);
                yMax = MAXVALUE(yMax, (float)myPoint.point->utm.y);
                valueMin = MINVALUE(valueMin, myPoint.value);
                valueMax = MAXVALUE(valueMax, myPoint.value);
            }

            interpolationPoints.push_back(myPoint);

            if (checkLapseRateCode(myPoint.lapseRateCode, mySettings->getUseLapseRateCode(), false))
                nrValid++;
        }
    }

    if (nrValid == 0)
        return false;

    mySettings->setPointsBoundingBoxArea((xMax - xMin) * (yMax - yMin));
    mySettings->setPointsRange(valueMin, valueMax);
    return true;
}
