#include "waterTableWidget.h"

WaterTableWidget::WaterTableWidget(const QString &id, std::vector<QDate> myDates, std::vector<float> myHindcastSeries,
                                   std::vector<float> myInterpolateSeries, QMap<QDate, float> obsDepths, float maxObservedDepth)
{
    this->setWindowTitle("Graph Id well: "+ id);
    this->resize(1240, 700);
    this->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    this->setAttribute(Qt::WA_DeleteOnClose);

    // layout
    QHBoxLayout *mainLayout = new QHBoxLayout();
    QVBoxLayout *plotLayout = new QVBoxLayout();

    waterTableChartView = new WaterTableChartView();
    plotLayout->addWidget(waterTableChartView);

    // menu
    QMenuBar* menuBar = new QMenuBar();
    QMenu *editMenu = new QMenu("Edit");

    menuBar->addMenu(editMenu);
    mainLayout->setMenuBar(menuBar);
    QAction* exportInterpolation = new QAction(tr("&Export interpolation as csv"), this);
    editMenu->addAction(exportInterpolation);

    mainLayout->addLayout(plotLayout);
    setLayout(mainLayout);

    connect(exportInterpolation, &QAction::triggered, this, &WaterTableWidget::on_actionExportInterpolationData);

    waterTableChartView->draw(myDates, myHindcastSeries, myInterpolateSeries, obsDepths, maxObservedDepth);
}

void WaterTableWidget::on_actionExportInterpolationData()
{

    QString csvFileName = QFileDialog::getSaveFileName(this, tr("Save current data"), "", tr("csv files (*.csv)"));
    if (csvFileName != "")
    {
        QFile myFile(csvFileName);
        if (!myFile.open(QIODevice::WriteOnly | QFile::Truncate))
        {
            QMessageBox::information(nullptr, "Error", "Open CSV failed: " + csvFileName + "\n ");
            return;
        }

        QTextStream myStream (&myFile);
        myStream.setRealNumberNotation(QTextStream::FixedNotation);
        myStream.setRealNumberPrecision(3);
        QString header = "date,value[m]";
        myStream << header << "\n";
        QDateTime firstDate(QDate(1970,1,1), QTime(0,0,0));
        QList<QPointF> dataPoints = waterTableChartView->exportInterpolationValues();
        for (int i = 0; i < dataPoints.size(); i++)
        {
            float xPointValue = dataPoints[i].x();
            QDateTime xValue = firstDate.addMSecs(xPointValue);
            QDate myDate = xValue.date().addDays(1);
            myStream << myDate.toString("yyyy-MM-dd") << "," << dataPoints[i].y()/100 << "\n";
        }

        myFile.close();

        return;
    }
}

WaterTableWidget::~WaterTableWidget()
{

}

void WaterTableWidget::closeEvent(QCloseEvent *event)
{
    event->accept();
}
