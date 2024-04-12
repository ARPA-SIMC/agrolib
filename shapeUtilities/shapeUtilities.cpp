#include "shapeUtilities.h"
#include "shapeHandler.h"
#include <QFile>
#include <QFileInfo>


/*! cloneShapeFile
 * \brief make a copy of shapefile (delete old version)
 * and return filename of the cloned shapefile
 */
QString cloneShapeFile(QString refFileName, QString newFileName)
{
    QFileInfo refFileInfo(refFileName);
    QFileInfo newFileInfo(newFileName);

    QString refFile = refFileInfo.absolutePath() + "/" + refFileInfo.baseName();
    QString newFile = newFileInfo.absolutePath() + "/" + newFileInfo.baseName();

    QFile::remove(newFile + ".dbf");
    QFile::copy(refFile +".dbf", newFile +".dbf");

    QFile::remove(newFile +".shp");
    QFile::copy(refFile +".shp", newFile +".shp");

    QFile::remove(newFile +".shx");
    QFile::copy(refFile +".shx", newFile +".shx");

    QFile::remove(newFile +".prj");
    QFile::copy(refFile +".prj", newFile +".prj");

    return newFile + ".shp";
}

/*! copyShapeFile
 * \brief make a copy of shapefile (keep original version)
 * and return filename of the cloned shapefile
 */
QString copyShapeFile(QString refFileName, QString newFileName)
{
    QFileInfo refFileInfo(refFileName);
    QFileInfo newFileInfo(newFileName);

    QString refFile = refFileInfo.absolutePath() + "/" + refFileInfo.baseName();
    QString newFile = newFileInfo.absolutePath() + "/" + newFileInfo.baseName();

    QFile::copy(refFile +".dbf", newFile +".dbf");
    QFile::copy(refFile +".shp", newFile +".shp");
    QFile::copy(refFile +".shx", newFile +".shx");
    QFile::copy(refFile +".prj", newFile +".prj");

    return(newFile + ".shp");
}


bool cleanShapeFile(Crit3DShapeHandler &shapeHandler)
{
    if (! shapeHandler.existRecordDeleted()) return true;

    QFileInfo fileInfo(QString::fromStdString(shapeHandler.getFilepath()));
    QString refFile = fileInfo.absolutePath() + "/" + fileInfo.baseName();
    QString tmpFile = refFile + "_temp";

    shapeHandler.packSHP(tmpFile.toStdString());
    shapeHandler.packDBF(tmpFile.toStdString());
    shapeHandler.close();

    QFile::remove(refFile + ".dbf");
    QFile::copy(tmpFile + ".dbf", refFile + ".dbf");
    QFile::remove(tmpFile + ".dbf");

    QFile::remove(refFile + ".shp");
    QFile::copy(tmpFile + ".shp", refFile + ".shp");
    QFile::remove(tmpFile + ".shp");

    QFile::remove(refFile + ".shx");
    QFile::copy(tmpFile + ".shx", refFile + ".shx");
    QFile::remove(tmpFile + ".shx");

    return shapeHandler.open(shapeHandler.getFilepath());
}


bool computeAnomaly(Crit3DShapeHandler *shapeAnomaly, Crit3DShapeHandler *shape1, Crit3DShapeHandler *shape2,
                    std::string id1, std::string id2, std::string field1, std::string field2, QString fileName, QString &errorStr)
{
    QString newShapeFileName = copyShapeFile(QString::fromStdString(shape1->getFilepath()), fileName);

    if (! shapeAnomaly->open(newShapeFileName.toStdString()))
    {
        errorStr = "Error in create/open new shapefile: " + newShapeFileName;
        return false;
    }

    // remove fields
    int nrFields = shape1->getFieldNumbers();
    for (int i = 0; i < nrFields; i++)
    {
        std::string fieldName = shape1->getFieldName(i);
        if (fieldName != id1)
        {
            int fieldPos = shapeAnomaly->getFieldPos(fieldName);
            if (fieldPos != -1)
            {
                if (! shapeAnomaly->removeField(fieldPos))
                {
                    errorStr = "Error in delete field: " + QString::fromStdString(fieldName);
                    return false;
                }
            }
        }
    }

    // add anomaly field
    if (! shapeAnomaly->addField("anomaly", FTDouble, 10, 1))
    {
        errorStr = "Error in create field 'anomaly'.";
        return false;
    }
    int anomalyPos = shapeAnomaly->getFieldPos("anomaly");

    // compute values
    int fieldPos1 = shape1->getFieldPos(field1);
    int fieldPos2 = shape2->getFieldPos(field2);

    for (int i = 0; i < shape1->getShapeCount(); i++)
    {
        double anomalyValue = shape2->getNumericValue(i, fieldPos2) - shape1->getNumericValue(i, fieldPos1);
        shapeAnomaly->writeDoubleAttribute(i, anomalyPos, anomalyValue);
    }

    return true;
}
