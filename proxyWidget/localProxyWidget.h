#ifndef LOCALPROXYWIDGET_H
#define LOCALPROXYWIDGET_H

#include <QtWidgets>
#include <QtCharts>
#include "chartView.h"
#include "meteoPoint.h"
#include "interpolationSettings.h"
#include "interpolationPoint.h"

class Crit3DLocalProxyWidget : public QWidget
{
    Q_OBJECT

public:
    Crit3DLocalProxyWidget(double x, double y, double zGrid, double zDEM, gis::Crit3DGisSettings gisSettings, Crit3DInterpolationSettings* interpolationSettings, Crit3DMeteoPoint* meteoPoints, int nrMeteoPoints, meteoVariable currentVariable, frequencyType currentFrequency, QDate currentDate, int currentHour, Crit3DQuality* quality,  Crit3DInterpolationSettings* SQinterpolationSettings, Crit3DMeteoSettings *meteoSettings, Crit3DClimateParameters *climateParam, bool checkSpatialQuality);
    ~Crit3DLocalProxyWidget();
    void closeEvent(QCloseEvent *event);
    void updateDateTime(QDate newDate, int newHour);
    void updateFrequency(frequencyType newFrequency);
    void changeProxyPos(const QString proxyName);
    void changeVar(const QString varName);
    void plot();
    void climatologicalLRClicked(int toggled);
    void modelLRClicked(int toggled);

private:
    double x;
    double y;
    double zGrid;
    double zDEM;
    gis::Crit3DGisSettings gisSettings;
    Crit3DInterpolationSettings* interpolationSettings;
    Crit3DQuality* quality;
    Crit3DInterpolationSettings* SQinterpolationSettings;
    Crit3DMeteoSettings *meteoSettings;
    Crit3DMeteoPoint* meteoPoints;
    Crit3DClimateParameters *climateParam;
    meteoVariable currentVariable;
    int nrMeteoPoints;
    bool checkSpatialQuality;
    frequencyType currentFrequency;
    QDate currentDate;
    int currentHour;
    std::vector <Crit3DInterpolationDataPoint> outInterpolationPoints;
    std::vector <Crit3DInterpolationDataPoint> subsetInterpolationPoints;
    QComboBox comboVariable;
    QComboBox comboAxisX;
    QCheckBox detrended;
    QCheckBox climatologicalLR;
    QCheckBox modelLR;
    QCheckBox stationWeights;
    QTextEdit r2;
    QTextEdit lapseRate;
    QTextEdit par0;
    QTextEdit par1;
    QTextEdit par2;
    QTextEdit par3;
    QTextEdit par4;
    QTextEdit par5;
    ChartView *chartView;
    meteoVariable myVar;
    int proxyPos;
    std::vector <QGraphicsTextItem*> weightLabels;

    Crit3DTime getCurrentTime();

signals:
    void closeLocalProxyWidget();
};


#endif // PROXYWIDGET_H
