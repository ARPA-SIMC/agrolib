#ifndef METEOWIDGET_H
#define METEOWIDGET_H

    #include <QtWidgets>
    #include <QtCharts>
    #include "meteo.h"
    #include "meteoPoint.h"
    #include "callout.h"

    class Crit3DMeteoWidget : public QWidget
    {
        Q_OBJECT

        public:
            Crit3DMeteoWidget(bool isGrid_, QString projectPath, Crit3DMeteoSettings* meteoSettings_);
            ~Crit3DMeteoWidget() override;

            int getMeteoWidgetID() const { return meteoWidgetID; }
            void setMeteoWidgetID(int value) { meteoWidgetID = value; }

            void setCurrentDate(QDate myDate) { currentDate = myDate; }

            void setIsEnsemble(bool value);
            bool getIsEnsemble() { return isEnsemble; }
            void setNrMembers(int value) { nrMembers = value; }

            void setDailyRange(QDate firstDate, QDate lastDate);
            void setHourlyRange(QDate firstDate, QDate lastDate);
            void setMonthlyRange(QDate firstDate, QDate lastDate);

            void addMeteoPointsEnsemble(Crit3DMeteoPoint mp);

            void drawMeteoPoint(Crit3DMeteoPoint mp, bool isAppend);
            void drawEnsemble();

    private:
            int meteoWidgetID;
            bool isGrid;
            bool isEnsemble;
            bool isInitialized;
            int nrMembers;

            QVector<Crit3DMeteoPoint> meteoPoints;
            QVector<Crit3DMeteoPoint> meteoPointsEnsemble;
            Crit3DMeteoSettings* meteoSettings;

            frequencyType currentFreq;
            QDate firstDailyDate;
            QDate lastDailyDate;
            QDate firstHourlyDate;
            QDate lastHourlyDate;
            QDate firstMonthlyDate;
            QDate lastMonthlyDate;
            QDate currentDate;

            QAction* dataSum;
            QPushButton *addVarButton;
            QPushButton *dailyButton;
            QPushButton *hourlyButton;
            QPushButton *monthlyButton;
            QPushButton *tableButton;
            QPushButton *redrawButton;
            QPushButton *shiftPreviousButton;
            QPushButton *shiftFollowingButton;
            QDateTimeEdit *firstDate;
            QDateTimeEdit *lastDate;
            QChartView *chartView;
            QChart *chart;
            QBarCategoryAxis *axisX;
            QBarCategoryAxis *axisXvirtual;
            QValueAxis *axisY_sx;
            QValueAxis *axisY_dx;
            QMap<QString, QList<QString>> MapCSVDefault;
            QMap<QString, QList<QString>> MapCSVStyles;
            QList<QString> currentVariables;
            QList<QString> nameLines;
            QList<QString> nameBar;
            double maxEnsembleBar;
            double maxEnsembleLine;
            double minEnsembleLine;
            QVector<QColor> colorLines;
            QMap<QString, QList<QColor>> colorLinesMpAppended;
            QVector<QColor> colorBar;
            QMap<QString, QList<QColor>> colorBarMpAppended;
            QLineSeries* zeroLine;
            QVector<QVector<QLineSeries*>> lineSeries;
            QVector<QBarSeries*> barSeries;
            QVector<QBoxPlotSeries*> ensembleSeries;
            QVector<QList<QBoxSet*>> ensembleSet;
            QVector<QVector<QBarSet*>> setVector;
            QList<QString> categories;
            QList<QString> categoriesVirtual;
            QList<QString> varToSumList;

            bool isLine;
            bool isBar;
            Callout *m_tooltip;

            void updateTimeRange();
            void resetValues();
            void resetEnsembleValues();
            void drawDailyVar();
            void drawEnsembleDailyVar();
            void drawHourlyVar();
            void drawMonthlyVar();
            void showMonthlyGraph();
            void showDailyGraph();
            void showHourlyGraph();
            void updateSeries();
            void redraw();
            void shiftPrevious();
            void shiftFollowing();
            void showTable();
            void showVar();
            void tooltipLineSeries(QPointF point, bool state);
            void editLineSeries();
            bool computeTooltipLineSeries(QLineSeries *series, QPointF point, bool state);
            void tooltipBar(bool state, int index, QBarSet *barset);
            void editBar();
            void handleMarkerClicked();
            void closeEvent(QCloseEvent *event) override;

            void on_actionChangeLeftAxis();
            void on_actionChangeRightAxis();
            void on_actionExportGraph();
            void on_actionRemoveStation();
            void on_actionInfoPoint();
            void on_actionDataAvailability();
            void on_actionDataSum();

            void drawAxisTitle();
            void drawSum();
            void checkExistingData();


    signals:
        void closeWidgetPoint(int);
        void closeWidgetGrid(int);

    };


    qreal findMedian(std::vector<double> sortedList, int begin, int end);


#endif // METEOWIDGET_H
