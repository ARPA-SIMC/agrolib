#ifndef WATERTABLECHARTVIEW_H
#define WATERTABLECHARTVIEW_H

    #include <QtWidgets>
    #include <QtCharts>
    #include "waterTable.h"
    #include "callout.h"

    class WaterTableChartView : public QChartView
    {
        public:
            WaterTableChartView(QWidget *parent = 0);

            void draw(std::vector<QDate> &myDates, std::vector<float> &myHindcastSeries, std::vector<float> &myInterpolateSeries,
                      QMap<QDate, float> obsDepths, float maximumObservedDepth);

            void tooltipObsDepthSeries(QPointF point, bool state);
            void tooltipLineSeries(QPointF point, bool state);
            void handleMarkerClicked();
            QList<QPointF> exportInterpolationValues();

        private:
            QScatterSeries* obsDepthSeries;
            QLineSeries* hindcastSeries;
            QLineSeries* interpolationSeries;

            QDateTimeAxis* axisX;
            QValueAxis* axisY;
            Callout *m_tooltip;
    };

#endif // WATERTABLECHARTVIEW_H
