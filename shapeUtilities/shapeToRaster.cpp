#include <float.h>
#include <math.h>

#include "shapeToRaster.h"
#include "commonConstants.h"
#include "gis.h"
#include "basicMath.h"


bool initializeRasterFromShape(const Crit3DShapeHandler &shapeHandler, gis::Crit3DRasterGrid &raster, double cellSize)
{
    int nrShape = shapeHandler.getShapeCount();
    if (nrShape <= 0)
    {
        // void shapefile
        return false;
    }

    gis::Crit3DRasterHeader header;
    ShapeObject object;
    Box<double> bounds;

    double ymin = DBL_MAX;
    double xmin = DBL_MAX;
    double ymax = DBL_MIN;
    double xmax = DBL_MIN;

    for (int i = 0; i < nrShape; i++)
    {
        shapeHandler.getShape(i, object);
        bounds = object.getBounds();
        ymin = MINVALUE(ymin, bounds.ymin);
        xmin = MINVALUE(xmin, bounds.xmin);
        ymax = MAXVALUE(ymax, bounds.ymax);
        xmax = MAXVALUE(xmax, bounds.xmax);
    }

    xmin = round(xmin * 10.) * 0.1;
    ymin = round(ymin * 10.) * 0.1;
    header.cellSize = cellSize;
    header.llCorner.x = xmin;
    header.llCorner.y = ymin;
    header.nrRows = int(floor((ymax - ymin) / cellSize)) + 1;
    header.nrCols = int(floor((xmax - xmin) / cellSize)) + 1;
    header.flag = NODATA;

    return raster.initializeGrid(header);
}


bool fillRasterWithShapeNumber(gis::Crit3DRasterGrid &raster, const Crit3DShapeHandler &shapeHandler)
{
    int nrShape = shapeHandler.getShapeCount();
    if (nrShape <= 0)
    {
        // void shapefile
        return false;
    }

    ShapeObject object;
    double x, y;
    Box<double> bounds;
    int r0, r1, c0, c1;

    raster.emptyGrid();

    for (int shapeIndex = 0; shapeIndex < nrShape; shapeIndex++)
    {
        shapeHandler.getShape(shapeIndex, object);

        // get bounds
        bounds = object.getBounds();
        gis::getRowColFromXY(*(raster.header), bounds.xmin, bounds.ymax, &r0, &c0);
        gis::getRowColFromXY(*(raster.header), bounds.xmax, bounds.ymin, &r1, &c1);

        // check bounds
        r0 = MAXVALUE(r0-1, 0);
        r1 = MINVALUE(r1+1, raster.header->nrRows -1);
        c0 = MAXVALUE(c0-1, 0);
        c1 = MINVALUE(c1+1, raster.header->nrCols -1);

        for (int row = r0; row <= r1; row++)
        {
            for (int col = c0; col <= c1; col++)
            {
                if (isEqual(raster.value[row][col], raster.header->flag))
                {
                    raster.getXY(row, col, x, y);
                    if (object.pointInPolygon(x, y))
                    {
                        raster.value[row][col] = float(shapeIndex);
                    }
                }
            }
        }
    }

    return true;
}


bool fillRasterWithField(gis::Crit3DRasterGrid &raster, Crit3DShapeHandler &shapeHandler, const std::string &fieldName)
{
    int nrShape = shapeHandler.getShapeCount();
    if (nrShape <= 0)
    {
        // void shapefile
        return false;
    }

    ShapeObject object;
    double x, y, fieldValue;
    Box<double> bounds;
    int r0, r1, c0, c1;

    int fieldIndex = shapeHandler.getDBFFieldIndex(fieldName.c_str());

    for (int shapeIndex = 0; shapeIndex < nrShape; shapeIndex++)
    {
        shapeHandler.getShape(shapeIndex, object);

        fieldValue = shapeHandler.getNumericValue(shapeIndex, fieldIndex);

        if (! isEqual(fieldValue, NODATA))
        {
            // get bounds
            bounds = object.getBounds();
            gis::getRowColFromXY(*(raster.header), bounds.xmin, bounds.ymax, &r0, &c0);
            gis::getRowColFromXY(*(raster.header), bounds.xmax, bounds.ymin, &r1, &c1);

            // check bounds
            r0 = MAXVALUE(r0-1, 0);
            r1 = MINVALUE(r1+1, raster.header->nrRows -1);
            c0 = MAXVALUE(c0-1, 0);
            c1 = MINVALUE(c1+1, raster.header->nrCols -1);

            for (int row = r0; row <= r1; row++)
            {
                for (int col = c0; col <= c1; col++)
                {
                    if (isEqual(raster.value[row][col], raster.header->flag))
                    {
                        raster.getXY(row, col, x, y);
                        if (object.pointInPolygon(x, y))
                        {
                            raster.value[row][col] = float(fieldValue);
                        }
                    }
                }
            }
        }
    }

    return true;
}


bool rasterizeShape(Crit3DShapeHandler &shapeHandler, gis::Crit3DRasterGrid &newRaster,
                    const std::string &field, double cellSize)
{
    if (! initializeRasterFromShape(shapeHandler, newRaster, cellSize))
        return false;

    if (field == "Shape ID")
    {
        if (! fillRasterWithShapeNumber(newRaster, shapeHandler))
            return false;
    }
    else
    {
        if (! fillRasterWithField(newRaster, shapeHandler, field))
            return false;
    }

    return true;
}


bool rasterizeShapeWithRef(const gis::Crit3DRasterGrid &refRaster, gis::Crit3DRasterGrid &newRaster,
                             Crit3DShapeHandler &shapeHandler, const std::string &fieldName)
{
    newRaster.initializeGrid(*(refRaster.header));

    int nrShape = shapeHandler.getShapeCount();
    if (nrShape <= 0)
    {
        // void shapefile
        return false;
    }

    int fieldIndex = NODATA;
    if (fieldName != "Shape ID")
    {
        fieldIndex = shapeHandler.getDBFFieldIndex(fieldName.c_str());
    }

    ShapeObject object;
    double x, y, fieldValue;
    Box<double> bounds;
    int r0, r1, c0, c1;

    for (int shapeIndex = 0; shapeIndex < nrShape; shapeIndex++)
    {
        shapeHandler.getShape(shapeIndex, object);

        if (fieldName == "Shape ID")
        {
            fieldValue = shapeIndex;
        }
        else
        {
            fieldValue = shapeHandler.getNumericValue(shapeIndex, fieldIndex);
        }

        if (! isEqual(fieldValue, NODATA))
        {
            // get bounds
            bounds = object.getBounds();
            gis::getRowColFromXY(*(newRaster.header), bounds.xmin, bounds.ymax, &r0, &c0);
            gis::getRowColFromXY(*(newRaster.header), bounds.xmax, bounds.ymin, &r1, &c1);

            // check bounds
            r0 = MAXVALUE(r0-1, 0);
            r1 = MINVALUE(r1+1, newRaster.header->nrRows -1);
            c0 = MAXVALUE(c0-1, 0);
            c1 = MINVALUE(c1+1, newRaster.header->nrCols -1);

            for (int row = r0; row <= r1; row++)
            {
                for (int col = c0; col <= c1; col++)
                {
                    if (! isEqual(refRaster.value[row][col], refRaster.header->flag))
                    {
                        newRaster.getXY(row, col, x, y);
                        if (object.pointInPolygon(x, y))
                        {
                            newRaster.value[row][col] = float(fieldValue);
                        }
                    }
                }
            }
        }
    }

    return true;
}

