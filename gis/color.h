#ifndef CRIT3DCOLOR_H
#define CRIT3DCOLOR_H

    #include <vector>

    namespace classificationMethod
    {
        enum type{EqualInterval, Gaussian, Quantile, Categories, UserDefinition };
    }

    class Crit3DColor {
    public:
        short red;
        short green;
        short blue;

        Crit3DColor();
        Crit3DColor(short, short, short);
    };

    class Crit3DColorScale {

    private:
        unsigned int _nrColors, _nrKeyColors;
        std::vector<Crit3DColor> color;
        double _minimum, _maximum;
        bool _isFixedRange;
        bool _isHideMinimum;
        bool _isTransparent;
        int _classification;

    public:
        std::vector<Crit3DColor> keyColor;

        Crit3DColorScale();

        void initialize(unsigned int nrKeyColors, unsigned int nrColors);
        bool classify();

        unsigned int nrColors() { return _nrColors; }
        unsigned int nrKeyColors() { return _nrKeyColors; }

        double minimum() { return _minimum; }
        void setMinimum(double min) { _minimum = min; }

        double maximum() { return _maximum; }
        void setMaximum(double max) { _maximum = max; }

        Crit3DColor* getColor(float myValue);
        unsigned int getColorIndex(float myValue);

        bool setRange(float minimum, float maximum);

        void setFixedRange(bool fixedRange) { _isFixedRange = fixedRange; }
        bool isFixedRange() { return _isFixedRange; }

        void setHideMinimum(bool isHideMinimum) { _isHideMinimum = isHideMinimum; }
        bool isHideMinimum() { return _isHideMinimum; }

        void setTransparent(bool transparent) { _isTransparent = transparent; }
        bool isTransparent() { return _isTransparent; }
    };

    bool setDefaultScale(Crit3DColorScale* myScale);
    bool setDTMScale(Crit3DColorScale* myScale);
    bool setTemperatureScale(Crit3DColorScale* myScale);
    bool setSlopeStabilityScale(Crit3DColorScale* myScale);
    bool setAnomalyScale(Crit3DColorScale* myScale);
    bool setPrecipitationScale(Crit3DColorScale* myScale);
    bool setRelativeHumidityScale(Crit3DColorScale* myScale);
    bool setRadiationScale(Crit3DColorScale* myScale);
    bool setWindIntensityScale(Crit3DColorScale* myScale);
    bool setCenteredScale(Crit3DColorScale* myScale);
    bool setCircolarScale(Crit3DColorScale* myScale);
    bool roundColorScale(Crit3DColorScale* myScale, int nrIntervals, bool lessRounded);
    bool reverseColorScale(Crit3DColorScale* myScale);
    bool setGrayScale(Crit3DColorScale* myScale);
    bool setBlackScale(Crit3DColorScale* myScale);
    bool setSurfaceWaterScale(Crit3DColorScale* myScale);
    bool setLAIScale(Crit3DColorScale* myScale);

    void mixColor(const Crit3DColor &backColor, const Crit3DColor &foreColor, Crit3DColor &colorOut, float alpha);

#endif // CRIT3DCOLOR_H
