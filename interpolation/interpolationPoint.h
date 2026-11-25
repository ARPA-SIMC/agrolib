#ifndef INTERPOLATIONPOINT_H
#define INTERPOLATIONPOINT_H

    #ifndef GIS_H
        #include "gis.h"
    #endif
    #ifndef METEO_H
        #include "meteo.h"
    #endif
    #ifndef INTERPOLATIONSETTINGS_H
        #include "interpolationSettings.h"
    #endif
    #include <memory>

    class Crit3DInterpolationDataPoint {
    private:

    public:
        std::shared_ptr<gis::Crit3DPoint> point;     // shared by the copies of the data point, freed with the last one

        int index;
        bool isActive;
        bool isMarked;
        float value;
        float regressionWeight;
        lapseRateCodeType lapseRateCode;
        gis::Crit3DRasterGrid* topographicDistanceMap;
        std::vector<float> proxyValues;

        Crit3DInterpolationDataPoint();

        float getProxyValue(unsigned int pos) const;
    };

#endif // INTERPOLATIONPOINT_H
