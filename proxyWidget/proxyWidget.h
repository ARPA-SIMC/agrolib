﻿#ifndef PROXYWIDGET_H
#define PROXYWIDGET_H

    #include <QtWidgets>
    #include <QtCharts>
    #include "chartView.h"
    #include "meteoPoint.h"
    #include "interpolationSettings.h"
    #include "interpolationPoint.h"

    class Crit3DProxyWidget : public QWidget
    {
        Q_OBJECT

        public:
        Crit3DProxyWidget(Crit3DInterpolationSettings* interpolationSettings, Crit3DMeteoPoint* meteoPoints, int nrMeteoPoints,
                          frequencyType currentFrequency, QDate currentDate, int currentHour, Crit3DQuality* quality,
                          Crit3DInterpolationSettings* SQinterpolationSettings, Crit3DMeteoSettings *meteoSettings,
                          Crit3DClimateParameters *climateParameters, bool checkSpatialQuality, int macroAreaNumber);

            ~Crit3DProxyWidget() override;
            void closeEvent(QCloseEvent *event) override;
            void updateDateTime(QDate newDate, int newHour);
            void updateFrequency(frequencyType newFrequency);
            void changeProxyPos(const QString proxyName);
            void changeVar(const QString varName);
            void plot();
            void climatologicalLRClicked(int toggled);
            void modelLRClicked(int toggled);
            void addMacroAreaLR();

    private:
            Crit3DInterpolationSettings* _interpolationSettings;
            Crit3DMeteoPoint* _meteoPoints;
            int _nrMeteoPoints;
            frequencyType _currentFrequency;
            QDate _currentDate;
            int _currentHour;
            Crit3DQuality* _quality;
            Crit3DInterpolationSettings* _SQinterpolationSettings;
            Crit3DMeteoSettings* _meteoSettings;

            bool _checkSpatialQuality;
            int _proxyPos;
            int _macroAreaNumber;

            std::vector <Crit3DInterpolationDataPoint> _outInterpolationPoints;
            Crit3DClimateParameters* _climateParameters;

            QValueAxis *axisY_sx;
            QComboBox comboVariable;
            QComboBox comboAxisX;
            QCheckBox detrended;
            QCheckBox climatologicalLR;
            QCheckBox modelLR;
            QTextEdit r2;
            QTextEdit lapseRate;
            ChartView *chartView;
            meteoVariable myVar;

            Crit3DTime getCurrentTime();

            void on_actionChangeLeftAxis();

    signals:
        void closeProxyWidget();
    };


#endif // PROXYWIDGET_H
