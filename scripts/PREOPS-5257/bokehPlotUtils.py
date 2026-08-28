import os

from astropy.table import vstack, Table
import math
import numpy as np
import pandas as pd

import bokeh
import bokeh.palettes as palettes
from bokeh.io import output_file, output_notebook, show, save, export_png
from bokeh.layouts import gridplot
from bokeh.models import (ColumnDataSource, CustomJS, Range1d, HoverTool,
                          Selection, TapTool, OpenURL, PrintfTickFormatter,
                          LinearColorMapper, ColorBar, NumeralTickFormatter)
from bokeh.models.annotations import Span, BoxAnnotation, Label
from bokeh.plotting import figure, reset_output
from bokeh.transform import factor_cmap, factor_mark, transform

from lsst.afw.cameraGeom import FOCAL_PLANE
from lsst.pipe.base import Struct

# This is a collection of plotting utilities, mostly for creating standalone
# interactive Bokeh plots of various visit-level, per-detector parameters and
# metrics.  It is currently being called from the makeVisitSummaryPlots.py
# script and a few notebooks in
# /sdf/data/rubin/user/laurenma/tickets/notebooks

# Set some global parameters.
N_DATA_HTML_MAX = 50000  # If nData is greater than this, only generate static png plots.
SORTED_FULL_BAND_LIST = ["u", "g", "r", "i", "z", "y", "N921"]
HOVER_COLOR = "cyan"
NON_SELECT_COLOR = "silver"


def setThreshParList():
    # Parameters with associated thresholds for inclusion in coadds.
    threshParList = [
        "astromOffsetMean",
        "medianE",
        "psfStarScaledDeltaSizeScatter",
        "psfTraceRadiusDelta",
        "psfTraceRadiusSigmaScaledDelta",
        "psfApFluxDelta",
        "psfApCorrSigmaScaledDelta",
        "psfSigma"
        # "psfFwhm"
    ]
    return threshParList


def setBandColorDict():
    bandColorDict = {
        "u": "#17becf",
        "g": "#1f77b4",
        "r": "#ff7f0e",
        "i": "#2ca02c",
        "z": "#d62728",
        "y": "#9467bd",
        "N921": "#8c564b",
        "N387": "#e377c2",
        "N816": "#bcbd22",
        # "": "#7f7f7f",
    }
    return bandColorDict


def setSimsBandColorDict():
    # This is the SIMS mapping.
    bandColorDict = {
        "u": "#4A7DB3",
        "g": "#68AD57",
        "r": "#EE8632",
        "i": "#E887BD",
        "z": "#D1352B",
        "y": "#8E529F",
        "N921": "#8c564b",
        "N387": "#e377c2",
        "N816": "#bcbd22",
    }
    return bandColorDict


def setDesBandColorDict():
    # This is the DES mapping...going with this one for now.
    bandColorDict = {
        "u": "#56B4E9",
        "g": "#008060",
        "r": "#FF4000",
        "i": "#850000",
        "z": "#6600CC",
        "y": "#000000",
        "N921": "#8c564b",
        "N387": "#e377c2",
        "N816": "#bcbd22",
    }
    return bandColorDict


def setBandMarkerDict():
    bandMarkerDict = {
        "u": "diamond",
        "g": "star",
        "r": "circle",
        "i": "square",
        "z": "hex",
        "y": "triangle",
        "N921": "plus",
        "N387": "square_pin",
        "N816": "triangle_pin",
    }
    return bandMarkerDict


def setUnitsDict():
    unitsDict=dict(
        astromOffsetMean="arcsec",
        astromOffsetStd="arcsec",
        psfSigma="pixels",
        apRad="pixels",
        psfFwhm="arcsec",
        seeing="arcsec",
        skyBg="counts",
        skyNoise="counts",
        visit="number",
        seqNum="number",
        detector="number",
        zenithDistance="degrees",
        zeroPoint="mag",
        expTimeScaledZp="mag",
        ra="degrees",
        dec="degrees",
        xFp="mm",
        yFp="mm",
        medianE="",
        psfStarScaledDeltaSizeScatter="",
        psfTraceRadiusDelta="pixels",
        psfTraceRadiusSigmaScaledDelta="pixels",
        psfApFluxDelta="normalized",
        psfApCorrSigmaScaledDelta="factor",
        maxDistToNearestPsf="pixles",
        effTime="sec",
        effTimePsfSigmaScale="",
        effTimeSkyBgScale="",
        effTimeZeroPointScale="",
        effTimeScale="t_eff = effTime/expTime",
        expTimeScaledSkyBg="counts/sec",
    )
    return unitsDict


def setAxisScaleDict():
    axisScaleDict = dict(
        astromOffsetMean="log",
        astromOffsetStd="log",
        psfSigma="linear",
        psfFwhm="linear",
        seeing="linear",
        apRad="linear",
        skyBg="linear",
        skyNoise="linear",
        visit="linear",
        seqNum ="linear",
        detector="linear",
        zenithDistance="linear",
        zeroPoint="linear",
        expTimeScaledZp="linear",
        ra="linear",
        dec="linear",
        xFp="linear",
        yFp="linear",
        medianE="linear",
        psfStarScaledDeltaSizeScatter="linear",
        psfTraceRadiusDelta="linear",
        psfTraceRadiusSigmaScaledDelta="linear",
        psfApFluxDelta="linear",
        psfApCorrSigmaScaledDelta="linear",
        maxDistToNearestPsf="linear",
        effTime="linear",
        effTimePsfSigmaScale="linear",
        effTimeSkyBgScale="linear",
        effTimeZeroPointScale="linear",
        effTimeScale="linear",
        expTimeScaledSkyBg="linear",
    )
    return axisScaleDict


def setAxisLimitsDict(instrument=None, threshStruct=None):
    axisLimitsDict = dict(
        astromOffsetMean=(0.006, 0.6),
        astromOffsetStd=(0.006, 0.6),
        psfSigma=(0.0, 5.5),
        psfFwhm=None,  # (0.3, 1.8),   #  None,  # (0.0, 1.99),
        seeing=None,
        apRad=None,
        skyBg=None,
        skyNoise=None,
        visit=None,
        detector=None,
        expTimeScaledSkyBg=None,
        zenithDistance=(0.0, 90.0),
        zeroPoint=(27.0, 34.0),
        expTimeScaledZp=(23.0, 28.0),
        medianE=(0.00002, 0.03),
        psfStarScaledDeltaSizeScatter=(0, 0.05),
        psfTraceRadiusDelta=(0.0, 8.0),
        psfTraceRadiusSigmaScaledDelta=(0.0, 3.0),
        psfApFluxDelta=(0.0, 2.0),
        psfApCorrSigmaScaledDelta=(0.0, 2.0),
        maxDistToNearestPsf=(0.0, 3000.0),
    )

    if instrument is not None:
        if instrument == "LSSTComCamSim":
            axisLimitsDict["psfSigma"] = (1.0, 12.0)
            # axisLimitsDict["psfFwhm"] = (0.0, 1.5)
            axisLimitsDict["zeroPoint"] = (27.0, 32.0)
            axisLimitsDict["expTimeScaledZp"] = (19.2, 24.8)
            axisLimitsDict["psfTraceRadiusDelta"] = (0.0, 8.0)
            axisLimitsDict["medianE"] = (0.00002, 0.12)
            axisLimitsDict["maxDistToNearestPsf"] = (0.0, 3200.0)

        if "LSSTCam" in instrument and "Sim" in instrument:
            axisLimitsDict["psfSigma"] = (0.0, 7.0)
            # axisLimitsDict["psfFwhm"] = (0.0, 2.2)
            axisLimitsDict["zeroPoint"] = (27.0, 32.0)
            axisLimitsDict["expTimeScaledZp"] = (26.5, 29.0)
            axisLimitsDict["psfTraceRadiusDelta"] = (0.0, 8.5)
            axisLimitsDict["medianE"] = (0.00002, 0.05)
            axisLimitsDict["maxDistToNearestPsf"] = (0.0, 3200.0)

        if instrument == "LSSTComCamSim":
            axisLimitsDict["psfSigma"] = (0.0, 3.8)
            # axisLimitsDict["psfFwhm"] = (0.0, 1.5)
            axisLimitsDict["expTimeScaledZp"] = (27.2, 27.8)
            axisLimitsDict["psfTraceRadiusDelta"] = (0.0, 0.9)
            axisLimitsDict["medianE"] = (0.00002, 0.008)

    if threshStruct is not None:
        fact = 5.0
        threshParList = setThreshParList()
        axisLimitsDict["astromOffsetMean"] = (axisLimitsDict["astromOffsetMean"][0],
                                              fact*threshStruct.maxMeanDistanceArcsec)
        axisLimitsDict["astromOffsetStd"] = axisLimitsDict["astromOffsetStd"]
        axisLimitsDict["medianE"] = (axisLimitsDict["medianE"][0],
                                     fact*threshStruct.maxEllipResidual)
        axisLimitsDict["psfStarScaledDeltaSizeScatter"] = (
            axisLimitsDict["psfStarScaledDeltaSizeScatter"][0],
            fact*threshStruct.maxScaledSizeScatter
        )
        axisLimitsDict["psfTraceRadiusDelta"] = (axisLimitsDict["psfTraceRadiusDelta"][0],
                                                 fact*threshStruct.maxPsfTraceRadiusDelta)
        axisLimitsDict["psfTraceRadiusDelta"] = (axisLimitsDict["psfTraceRadiusSigmaScaledDelta"][0],
                                                 fact*threshStruct.maxPsfTraceRadiusSigmaScaledDelta)
        axisLimitsDict["psfApFluxDelta"] = (axisLimitsDict["psfApFluxDelta"][0],
                                            fact*threshStruct.maxPsfApFluxDelta)
        axisLimitsDict["psfApCorrSigmaScaledDelta"] = (axisLimitsDict["psfApCorrSigmaScaledDelta"][0],
                                                       fact*threshStruct.maxPsfApCorrSigmaScaledDelta)
    return axisLimitsDict


def introspectConfigs(instrument=None, butler=None):
    # Stack Defaults:
    maxMeanDistanceArcsec = 0.5
    maxEllipResidual = 0.007
    maxScaledSizeScatter = 0.019
    maxPsfTraceRadiusDelta = 0.7
    maxPsfTraceRadiusSigmaScaledDelta = 0.3
    maxPsfApFluxDelta = 0.24
    maxPsfApCorrSigmaScaledDelta = 0.22

    if butler is not None:
        try:
            # Introspect config value from persisted configs.
            calibrateConfigDatarefList = list(butler.registry.queryDatasets("calibrate_config"))
            calibConfigFile = butler.get(calibrateConfigDatarefList[0])
            maxMeanDistanceArcsec = calibConfigFile.astrometry.maxMeanDistanceArcsec
        except LookupError:
            print("Could not load calibrate_config...setting fiducial or known override for "
                  "maxMeanDistanceArcsec.")
            if instrument == "LSSTCam-imSim":
                maxMeanDistanceArcsec = 0.05

        try:
            makeWarpConfigDatarefList = list(butler.registry.queryDatasets("makeWarp_config"))
            makeWarpConfigFile = butler.get(makeWarpConfigDatarefList[0])
            maxEllipResidual = makeWarpConfigFile.select.maxEllipResidual
            maxScaledSizeScatter = makeWarpConfigFile.select.maxScaledSizeScatter
            maxPsfTraceRadiusDelta = makeWarpConfigFile.select.maxPsfTraceRadiusDelta
        except LookupError:
            print("Could not load makeWarp_config...setting fiducial or known overrides for "
                  "maxEllipResidual, maxScaledSizeScatter, and maxPsfTraceRadiusDelta.")
            if instrument == "LSSTCam-imSim":
                maxMeanDistanceArcsec = 0.05
                maxEllipResidual = 0.0036
                maxScaledSizeScatter = 0.011
                maxPsfTraceRadiusDelta = 0.09
                maxPsfApFluxDelta = 0.062
                maxPsfApCorrSigmaScaledDelta = 0.053
            if instrument == "LSSTComCamSim":
                maxEllipResidual = 0.0026
                maxScaledSizeScatter = 0.014
                maxPsfTraceRadiusDelta = 0.065
                maxPsfApFluxDelta = 0.047
                maxPsfApCorrSigmaScaledDelta = 0.041
            elif instrument == "LATISS":
                maxEllipResidual = 0.027
                maxScaledSizeScatter = 0.026
                maxPsfTraceRadiusDelta = 2.9
                maxPsfApFluxDelta = 0.075
                maxPsfApCorrSigmaScaledDelta = 0.118
    else:
        print("No butler provided, so can't introspect configs.  Setting fiducal or known overrides "
              "for all thresholds.")
        if instrument == "LSSTCam-imSim":
            maxMeanDistanceArcsec = 0.05
            maxEllipResidual = 0.0036
            maxScaledSizeScatter = 0.011
            maxPsfTraceRadiusDelta = 0.09
            maxPsfApFluxDelta = 0.062
            maxPsfApCorrSigmaScaledDelta = 0.053
        elif instrument == "LSSTComCamSim":
            maxEllipResidual = 0.0026  # 0.004
            maxScaledSizeScatter = 0.014
            maxPsfTraceRadiusDelta = 0.065  # 0.091
            maxPsfApFluxDelta = 0.047
            maxPsfApCorrSigmaScaledDelta = 0.041
        elif instrument == "LATISS":
            maxEllipResidual = 0.027
            maxScaledSizeScatter = 0.026
            maxPsfTraceRadiusDelta = 2.9
            maxPsfTraceRadiusSigmaScaledDelta = 1.1
            maxPsfApFluxDelta = 0.075
            maxPsfApCorrSigmaScaledDelta = 0.118

    return Struct(
        maxMeanDistanceArcsec=maxMeanDistanceArcsec,
        maxEllipResidual=maxEllipResidual,
        maxScaledSizeScatter=maxScaledSizeScatter,
        maxPsfTraceRadiusDelta=maxPsfTraceRadiusDelta,
        maxPsfTraceRadiusSigmaScaledDelta=maxPsfTraceRadiusSigmaScaledDelta,
        maxPsfApFluxDelta=maxPsfApFluxDelta,
        maxPsfApCorrSigmaScaledDelta=maxPsfApCorrSigmaScaledDelta,
    )


def loadVisitData(butler, camera):
    # Figure out if ccdVisitTable datasets exist.  If not, load visitSummary tables.
    doLoadCcdVisitTables, doLoadVisitSummaries = False, False
    if doTryCcdVisitTable:
        try:
            butler.exists("ccdVisitTable")
            doLoadCcdVisitTables = True
            dataSourceStr = "ccdVisitTable"
            print("Found ccdVisitTable dataset, so will load those...")
        except LookupError:
            doLoadVisitSummaries = True
            dataSourceStr = "visitSummary"
            print("Could not find ccdVisitTable dataset, so will load visitSummay tables...")
    else:
        doLoadVisitSummaries = True
        dataSourceStr = "visitSummary"
    if doLoadVisitSummaries:
        visitDataAll = loadVisitSummaries(butler)
    elif doLoadCcdVisitTables:
        visitDataAll = loadCcdVisitTables(butler)
    else:
        print("No visit table datasets found (looked for ccdVisitTable & visitSummay)")
    print("dataSourceStr = {}".format(dataSourceStr))
    visitDataAll = addDetCenter(camera, visitDataAll)

    return visitDataAll, dataSourceStr


def loadCcdVisitTables(butler, visitVetoList=[]):
    # Read in data from ccdVisitTable when available
    ccdVisitTable = butler.get("ccdVisitTable")
    for visit in visitVetoList:
        if visit in ccdVisitTable.visitId.values:
            ccdVisitTable = ccdVisitTable[ccdVisitTable.visitId != visit]
    ccdVisitTable["medianE"] = np.sqrt(ccdVisitTable["psfStarDeltaE1Median"]**2.
                                       + ccdVisitTable["psfStarDeltaE2Median"]**2.)
    ccdVisitTable["visit"] = ccdVisitTable["visitId"]
    ccdVisitTable["expTimeScaledZp"] = ccdVisitTable["zeroPoint"] - 2.5*np.log10(ccdVisitTable["expTime"])
    return ccdVisitTable


def loadVisitSummaries(butler, whereStr=""):
    # Read in visitSummary tables instead of ccdVisitTable (e.g. if not available)
    # Load and stack all visitSummary tables in repo collection
    # print("visitSummary dataId search string:\n{}\n".format(whereStr))
    visitSummaryRefList = list(butler.registry.queryDatasets("visitSummary", where=whereStr))
    nVisitCollection = len(visitSummaryRefList)
    collectionsStr = ""
    for collect in butler.collections:
        collectionsStr += collect + " "
    print("Number of visits in collection {}: {}".format(collectionsStr, nVisitCollection))

    visitSummaryList = []
    nDigitVisit = len(str(nVisitCollection))
    maxNumVisit = nVisitCollection

    for iv, visitSummaryRef in enumerate(visitSummaryRefList[:maxNumVisit]):
        visitSummary = butler.get(visitSummaryRef)
        if "expTime" not in visitSummary.columns:
            expTime = visitSummary[0].getVisitInfo().exposureTime
        visitSummary = visitSummary.asAstropy()
        visitSummaryList.append(visitSummary)
        visit = visitSummaryRef.dataId["visit"]
        if "expTime" not in visitSummary.columns:
            visitSummary["expTime"] = [expTime]*len(visitSummary)

    visitSummaryAll = vstack(visitSummaryList)
    visitSummaryAll["detector"] = visitSummaryAll["id"]
    havePsfStats = False
    if "psfStarDeltaE1Median" in visitSummaryAll.columns:
        print("visitSummary does contain PSF metrics")
        havePsfStats = True
    if havePsfStats:
        visitSummaryAll["medianE"] = np.sqrt(visitSummaryAll["psfStarDeltaE1Median"]**2.
                                             + visitSummaryAll["psfStarDeltaE2Median"]**2.)
    visitSummaryAll["expTimeScaledZp"] = (visitSummaryAll["zeroPoint"]
                                          - 2.5*np.log10(visitSummaryAll["expTime"]))
    return visitSummaryAll


def addDetCenter(camera, visitDataAll):
    # Add the detector center in Focal Plane coords to visitDataAll.
    nDetectorCamera = len(camera)
    centerFpList = [camera[det].getCenter(FOCAL_PLANE) for det in visitDataAll["detector"]]
    xFpList = [centerFp[0] for centerFp in centerFpList]
    yFpList = [centerFp[1] for centerFp in centerFpList]
    visitDataAll["xFp"] = xFpList
    visitDataAll["yFp"] = yFpList
    return visitDataAll


def getMetricThresh(xCol, threshParList, threshStruct):
    xThresh = None
    if xCol in threshParList:
        if xCol == "astromOffsetMean":
            xThresh = threshStruct.maxMeanDistanceArcsec
        elif xCol == "medianE":
            xThresh = threshStruct.maxEllipResidual
        elif xCol == "psfStarScaledDeltaSizeScatter":
            xThresh = threshStruct.maxScaledSizeScatter
        elif xCol == "psfTraceRadiusDelta":
            xThresh = threshStruct.maxPsfTraceRadiusDelta
        elif xCol == "psfApFluxDelta":
            xThresh = threshStruct.maxPsfApFluxDelta
        elif xCol == "psfApCorrSigmaScaledDelta":
            xThresh = threshStruct.maxPsfApCorrSigmaScaledDelta
        else:
            print("Don't have threshold for {}".format(xCol))
    return xThresh


def getSource(dataTable, xColList, yColList, zColList=None,
              detectorStr="detector", bandList=None):
    # Function to create a column data source for the plots to share.
    # Sort by order in bandList for legend consistency.
    if bandList is not None:
        if isinstance(dataTable, Table):
            dataTable.remove_columns(["raCorners", "decCorners"])
            dataTable = dataTable.to_pandas()
        bandOrder = pd.api.types.CategoricalDtype(bandList, ordered=True)
        dataTable["band"] = dataTable["band"].astype(bandOrder)
        dataTable = dataTable.sort_values("band")

    dataDict = dict(visit=dataTable["visit"],
                    detector=dataTable[detectorStr],
                    band=dataTable["band"])
    for xCol, yCol in zip(xColList, yColList):
        dataDict[xCol] = dataTable[xCol]
        dataDict[yCol] = dataTable[yCol]
    if zColList:
        for zCol in zColList:
            dataDict[zCol] = dataTable[zCol]
    source = ColumnDataSource(data=dataDict)
    return source


# Create a custom hover tool on both panels.
def createBokehTools(TOOLTIPS_LIST, MENU_TOOLS):
    HOVER_TOOL = HoverTool(tooltips=TOOLTIPS_LIST)
    BOKEH_TOOLS = [MENU_TOOLS, HOVER_TOOL]  # custom_hover
    return BOKEH_TOOLS


def addBokehThresholds(xCol, yCol, bokehPlot, threshStruct, xData=None, yData=None,
                       expTime=30.0, plotThreshLine=False):
    xThresh = None
    yThresh = None
    if xCol == "astromOffsetMean":
        xThresh = threshStruct.maxMeanDistanceArcsec
    if yCol == "astromOffsetMean":
        yThresh = threshStruct.maxMeanDistanceArcsec
    if xCol == "medianE":
        xThresh = threshStruct.maxEllipResidual
    if yCol == "medianE":
        yThresh = threshStruct.maxEllipResidual
    if xCol == "psfStarScaledDeltaSizeScatter":
        xThresh = threshStruct.maxScaledSizeScatter
    if yCol == "psfStarScaledDeltaSizeScatter":
        yThresh = threshStruct.maxScaledSizeScatter
    if xCol == "psfTraceRadiusDelta":
        xThresh = threshStruct.maxPsfTraceRadiusDelta
    if yCol == "psfTraceRadiusDelta":
        yThresh = threshStruct.maxPsfTraceRadiusDelta
    if xCol == "psfApFluxDelta":
        xThresh = threshStruct.maxPsfApFluxDelta
    if yCol == "psfApFluxDelta":
        yThresh = threshStruct.maxPsfApFluxDelta
    if xCol == "psfApCorrSigmaScaledDelta":
        xThresh = threshStruct.maxPsfApCorrSigmaScaledDelta
    if yCol == "psfApCorrSigmaScaledDelta":
        yThresh = threshStruct.maxPsfApCorrSigmaScaledDelta
    if xCol == "psfApFluxDelta":
        xThresh = threshStruct.maxPsfApFluxDelta
    if yCol == "psfApFluxDelta":
        yThresh = threshStruct.maxPsfApFluxDelta

    if plotThreshLine:
        threshKwargs = dict(line_color="dimgray", line_width=2, line_dash=(10, 5))
    else:
        threshKwargs = dict(line_color=None)
    if yCol.startswith("effTime"):
        if yCol == "effTime":
            hLine = Span(location=expTime, dimension="width", **threshKwargs)
        else:
            hLine = Span(location=1.0, dimension="width", **threshKwargs)
        bokehPlot.add_layout(hLine)

    if xThresh is not None:
        xSpanThresh = Span(dimension="height", location=xThresh, **threshKwargs)
        bokehPlot.add_layout(xSpanThresh)
        threshBox = BoxAnnotation(left=xThresh, fill_alpha=0.1, fill_color="black")
        bokehPlot.add_layout(threshBox)
    if yThresh is not None:
        ySpanThresh = Span(dimension="width", location=yThresh, **threshKwargs)
        bokehPlot.add_layout(ySpanThresh)
        threshBox = BoxAnnotation(bottom=yThresh, fill_alpha=0.1, fill_color="black")
        bokehPlot.add_layout(threshBox)

    xNumRej = None
    yNumRej = None
    totNumRej = 0
    if xData is not None and xThresh is not None:
        xNumRej = sum(xData > xThresh) + sum(~np.isfinite(xData))
    if yData is not None and yThresh is not None:
        yNumRej = sum(yData > yThresh) + sum(~np.isfinite(yData))
    if xNumRej is not None and yNumRej is not None:
        for xDatum, yDatum in zip(xData, yData):
            if xDatum > xThresh or yDatum > yThresh:
                totNumRej += 1
    elif xNumRej is not None:
        totNumRej = xNumRej
    elif yNumRej is not None:
        totNumRej = yNumRej

    return bokehPlot, xNumRej, yNumRej, totNumRej


def addMetricColumns(dataTable, expTime=30.0, pixelScale=0.2):
    dataTable["psfTraceRadiusSigmaScaled"] = dataTable["psfTraceRadiusDelta"]/dataTable["psfSigma"]
    dataTable["medianE"] = np.sqrt(dataTable["psfStarDeltaE1Median"]**2.
                                + dataTable["psfStarDeltaE2Median"]**2.)
    if "expTime" not in dataTable.columns:
        print("Adding nominal expTime of {}s to dataTable".format(expTime))
        dataTable["expTime"] = expTime  # TODO: get this into the summaryStats!!
    dataTable["expTimeScaledZp"] = dataTable["zeroPoint"] - 2.5*np.log10(dataTable["expTime"])
    dataTable["expTimeScaledSkyBg"] = dataTable["skyBg"]/dataTable["expTime"]
    if "effTime" not in dataTable.columns:
        print("No effTime metrics...adding as NaNs")
        dataTable["effTime"] = float("nan")
        dataTable["effTimeScale"] = float("nan")
        dataTable["effTimePsfSigmaScale"] = float("nan")
        dataTable["effTimeSkyBgScale"] = float("nan")
        dataTable["effTimeZeroPointScale"] = float("nan")
        dataTable["expTimeScaledSkyBg"] = float("nan")
    else:
        dataTable["effTimeScale"] = dataTable["effTime"]/expTime

    if "pixelScale" not in dataTable.columns:
        print("Adding nominal pixelScale of {} arcsec/pixel to dataTable".format(pixelScale))
        dataTable["pixelScale"] = pixelScale  # TODO: get this into the summaryStats!!
    dataTable["psfFwhm"] = np.sqrt(8.0*np.log(2.0))*dataTable["psfSigma"]*dataTable["pixelScale"]
    return dataTable


def addSeqNumColumn(dataTable, visitOffset=None):
    # In the context of per-night plots, this should be the seq_num of the
    # exposure.  For other contexts, a visitOffset (integer) may make more
    # sense.  The logic here assumes the dataTable is sorted by visit,
    # detector.
    visitListInt = [int(visit) for visit in dataTable["visit"]]
    if visitOffset is None:
        visitNumListInt = [int(visit[8:]) for visit in dataTable["visit"]]
        maxVisitNum = max(visitNumListInt)
        dayNumPlus = 1000 if maxVisitNum < 1000 else 10000
        visitOffsetList = []
        baseDateVisit = min(dataTable["visit"])
        baseDateStr = baseDateVisit[:8]
        baseVisitOffset = int(baseDateStr + "00000")
        for visit in dataTable["visit"]:
            visitDateStr = visit[:8]
            currentVisitOffset = int(visitDateStr + "00000")
            deltaDay = int(visitDateStr) - int(baseDateStr)
            if deltaDay > 0:
                visitOffsetList.append(baseVisitOffset - deltaDay*dayNumPlus + deltaDay*100000)
            else:
                visitOffsetList.append(baseVisitOffset)
    else:
        visitOffsetList = len(visitListInt)*[visitOffset]
    seqNumList = [
        visitInt - visitOffset for (visitInt, visitOffset) in zip(visitListInt, visitOffsetList)
        ]
    dataTable["seqNum"] = seqNumList
    return dataTable


def makeParameterCompBokehPlots(visitDataAll, fileExtension, xColList, yColList,
                                plotTitleStr, plotNameStr, dataSourceStr,
                                saveDir, bandList, threshStruct, instrument,
                                collection=None, setXRange=False, setYRange=False):
    # Create figure of various parameter vs. parameter scatter plots.
    # If number of data points is not too large (>N_DATA_HTML_MAX) create
    # a linked 3x3 grid (doAddTools == True), otherwise persist as a
    # static png.

    # Set some parameter defaults for Bokeh plotting
    nBins = 50
    plotWidth = 520  # 350
    axisScaleDict = setAxisScaleDict()
    bandColorDict = setDesBandColorDict()
    palette = [bandColorDict[band] for band in bandList]
    index_cmap = factor_cmap("band", palette=palette, factors=bandList)
    bandMarkerDict = setBandMarkerDict()
    bandMarkers = [bandMarkerDict[band] for band in bandList]
    index_mark = factor_mark("band", bandMarkers, bandList)
    if setXRange or setYRange:
        axisLimitsDict = setAxisLimitsDict(instrument=instrument, threshStruct=threshStruct)

    nDataId = len(visitDataAll)
    visitList = list(set(visitDataAll["visit"]))
    if len(str(visitList[0])) == 13:
        obsDateFirst = str(min(visitList))[:8]
        obsDateLast = str(max(visitList))[:8]
    else:
        obsDateFirst = None
        obsDateLast = None
    nVisit = len(visitList)

    if fileExtension == "html":
        print("Creating interactive standalone html plot (nData = {}).".format(nDataId))
    else:
        print("Creating static png plot (nData = {})".format(nDataId))

    dataTable = visitDataAll.copy()
    sourceAll = getSource(dataTable, xColList, yColList, bandList=bandList)

    if fileExtension == "html":
        callback = CustomJS(
            args=dict(source=sourceAll),
            code="""
            var inds = source.selected.indices;
            if (inds.length > 1) {
                var alertStr = "Non-unique tap selection (N = " + inds.length + ").  Zoom in to select single points."
            }
            else {
                var alertStr = "visit: " + source.data.visit[inds] + "  detector: " + source.data.detector[inds]+ "  band: " + source.data.band[inds];
            }
            alert(alertStr);
            """)
        tapTool = TapTool(callback=callback)
        TOOLTIPS_LIST = [("visit", "@visit"), ("detector", "@detector"), ("band", "@band")]
        MENU_TOOLS = "box_select, lasso_select, box_zoom, reset, help"
        BOKEH_TOOLS = createBokehTools(TOOLTIPS_LIST, MENU_TOOLS)
    else:
        MENU_TOOLS, BOKEH_TOOLS, tapTool = None, None, None

    gridList = []
    iCol = -1

    for xCol, yCol in zip(xColList, yColList):
        iCol += 1
        xAxScale = axisScaleDict[xCol]
        yAxScale = axisScaleDict[yCol]
        titleStr = ""
        xRange, yRange = None, None
        if setXRange:
            if xCol in axisLimitsDict:
                xRange = axisLimitsDict[xCol]
            # else:
            #    print("{} not found in axisLimitsDict.  Not setting x-axis limits...".format(xCol))
        if setYRange:
            if yCol in axisLimitsDict:
                yRange = axisLimitsDict[yCol]
            # else:
            #     print("{} not found in axisLimitsDict.  Not setting y-axis limits...".format(yCol))
        if iCol == 0:
            titleStr = "{}".format(plotTitleStr)
        if iCol == 1:
            titleStr = "N_visit = {}   N_dataId = {}".format(nVisit, nDataId)
        if iCol == 2:
            if obsDateFirst is not None:
                titleStr = "dayObsStart: {}  dayObsEnd: {}".format(obsDateFirst, obsDateLast)
            elif collection is not None:
                titleStr = "{}".format(collection)
        if fileExtension == "html":
            if xRange is not None and yRange is not None:
                p = figure(tools=[MENU_TOOLS], width=plotWidth, height=plotWidth,
                           x_axis_type=xAxScale, y_axis_type=yAxScale,
                           title=titleStr, x_range=xRange, y_range=yRange)
            elif xRange is not None and yRange is None:
                p = figure(tools=[MENU_TOOLS], width=plotWidth, height=plotWidth,
                           x_axis_type=xAxScale, y_axis_type=yAxScale,
                           title=titleStr, x_range=xRange)
            elif xRange is None and yRange is not None:
                p = figure(tools=[MENU_TOOLS], width=plotWidth, height=plotWidth,
                           x_axis_type=xAxScale, y_axis_type=yAxScale,
                           title=titleStr, y_range=yRange)
            else:
                p = figure(tools=[MENU_TOOLS], width=plotWidth, height=plotWidth,
                           x_axis_type=xAxScale, y_axis_type=yAxScale,
                           title=titleStr)
            p.tools.append(tapTool)
        else:
            if xRange is not None and yRange is not None:
                p = figure(width=plotWidth, height=plotWidth,
                           x_axis_type=xAxScale, y_axis_type=yAxScale,
                           title=titleStr, x_range=xRange, y_range=yRange)
            elif xRange is not None and yRange is None:
                p = figure(width=plotWidth, height=plotWidth,
                           x_axis_type=xAxScale, y_axis_type=yAxScale,
                           title=titleStr, x_range=xRange)
            elif xRange is None and yRange is not None:
                p = figure(width=plotWidth, height=plotWidth,
                           x_axis_type=xAxScale, y_axis_type=yAxScale,
                           title=titleStr, y_range=yRange)
            else:
                p = figure(width=plotWidth, height=plotWidth,
                           x_axis_type=xAxScale, y_axis_type=yAxScale,
                           title=titleStr)
            p.toolbar.logo = None
            p.toolbar_location = None

        p.scatter(xCol, yCol, source=sourceAll, alpha=1.0, size=6, legend_field="band",
                  color=index_cmap,
                  marker=index_mark,
                  nonselection_fill_color=NON_SELECT_COLOR,
                  nonselection_line_color=NON_SELECT_COLOR)

        p, xNumRej, yNumRej, totNumRej = addBokehThresholds(
            xCol, yCol, p, threshStruct, xData=dataTable[xCol], yData=dataTable[yCol]
        )
        xStart = p.x_range.start
        yStart = p.y_range.start
        xDelta = p.x_range.end - xStart
        yDelta = p.y_range.end - yStart
        if not np.isfinite(xDelta):
            xData = dataTable[xCol]  # .copy(deep=True)
            if xCol == "detector":
                xData = pd.to_numeric(xData, errors="coerce")
            xDelta = 1.05*(np.nanmax(xData) - np.nanmin(xData))
            xStart = np.nanmin(xData)
        if not np.isfinite(yDelta):
            yDelta = 1.05*(np.nanmax(dataTable[yCol]) - np.nanmin(dataTable[yCol]))
            yStart = np.nanmin(dataTable[yCol])

        rejColor = "royalblue"
        nRejKwargs = dict(x_units="data", y_units="data",
                          text_font_size="15px", text_color="white",
                          border_line_color=rejColor, border_line_alpha=1.0,
                          background_fill_color=rejColor, background_fill_alpha=1.0,
                          border_line_width=8,
                          )
        if xNumRej is not None and yNumRej is not None:
            xNumRejLabel = Label(
                x=xStart + 0.98*xDelta, y=yStart + 0.05*yDelta,
                text="xNumRej = {}".format(xNumRej), text_align="right", **nRejKwargs
            )
            p.add_layout(xNumRejLabel)
            yNumRejLabel = Label(
                x=xStart + 0.02*xDelta, y=yStart + yDelta - len(bandList)*0.08*yDelta,
                text="yNumRej = {}".format(yNumRej), text_align="left", **nRejKwargs
            )
            p.add_layout(yNumRejLabel)
        if totNumRej > 0:
            totNumRejLabel = Label(
                x=xStart + 0.98*xDelta, y=yStart + 0.92*yDelta,
                text="totNumRej = {}".format(totNumRej), text_align="right", **nRejKwargs
            )
            p.add_layout(totNumRejLabel)

        if len(gridList) == 0:
            p.legend.location = "top_left"
            p.legend.click_policy = "hide"
            p.legend.padding = 0
            p.legend.spacing = 0
        else:
            p.legend.visible = False

        unitsDict = setUnitsDict()
        xUnitFormat = "%.2f"
        if unitsDict[xCol] in ["number", "degrees", "counts"]:
            xUnitFormat = "%.1e" if xCol == "visit" else "%d"
        if xCol in ["psfStarScaledDeltaSizeScatter", "medianE"] or "astrom" in xCol:
            xUnitFormat = "%.3f"
        yUnitFormat = "%.2f"
        if unitsDict[yCol] in ["number", "degrees", "counts"]:
            yUnitFormat = "%.1e" if yCol == "visit" else "%d"
        if yCol in ["psfStarScaledDeltaSizeScatter", "medianE"] or "astrom" in yCol:
            yUnitFormat = "%.3f"


        # if unitsDict[xCol] == "":
        #     unitFormat = "%.4f"
        # if unitsDict[xCol] in ["number", "degrees", "counts"]
        #     unitFormat = "%.1e" if xCol == "visit" else "%d"
        p.xaxis[0].formatter = PrintfTickFormatter(format=xUnitFormat)
        # unitFormat = "%.2f"
        # if unitsDict[yCol] == "":
        #     unitFormat = "%.4f"
        # if unitsDict[yCol] in ["number", "degrees", "counts"]:
        #    unitFormat = "%.1e" if yCol == "visit" else "%d"
        if yCol.startswith("effTime"):
            yDataMax = max(dataTable[yCol])
            yDataMin = min(dataTable[yCol])
            deltaLim = 2.0 if yCol == "effTime" else 0.05
            yAxisMax = max(1.2, yDataMax + deltaLim)
            yAxisMin = min(0.8, max(0.0, yDataMin - deltaLim))
            # if yDataMax >= 1.0:
            #     yAxisMin = min(0.8, max(0.0, 2.0 - yAxisMax))
            # else:
            #     yAxisMin = 0.0
            p.y_range = Range1d(yAxisMin, yAxisMax)
        p.yaxis[0].formatter = PrintfTickFormatter(format=yUnitFormat)
        p.xaxis.major_label_orientation = np.pi/4
        if unitsDict[xCol] == "":
            p.xaxis.axis_label = "{}".format(xCol, unitsDict[xCol])
        else:
            p.xaxis.axis_label = "{} ({})".format(xCol, unitsDict[xCol])
        if unitsDict[yCol] == "":
            p.yaxis.axis_label = "{}".format(yCol, unitsDict[yCol])
        else:
            p.yaxis.axis_label = "{} ({})".format(yCol, unitsDict[yCol])

        p.title.text_font_size = "13pt"
        p.xaxis.axis_label_text_font_size = "11pt"
        p.yaxis.axis_label_text_font_size = "11pt"

        if iCol == 0:
            gridList.append([p])
        elif np.mod(iCol, 3) > 0:
            gridList[-1].extend([p])
        else:
            gridList.append([p])

    fullPlot = gridplot(gridList)

    subSaveDir = "correlationPlots/"
    if not os.path.exists(saveDir + subSaveDir):
        print("Creating new directory: {}".format(saveDir + subSaveDir))
        os.mkdir(saveDir + subSaveDir)
    plotFilename = "{}{}{}_{}_{}.{}".format(
        saveDir, subSaveDir, instrument, dataSourceStr, plotNameStr, fileExtension)
    if fileExtension == "html":
        output_file(plotFilename,
                    title="{} {}: {}".format(instrument, dataSourceStr, plotTitleStr),
                    mode="inline")
        print("Saving {}".format(plotFilename))
        save(fullPlot)
    else:
        fullPlot.toolbar_location = None
        export_png(fullPlot, filename=plotFilename)
        print("Saving {}".format(plotFilename))
    return None


def makeScatAndHistPlots(visitDataAll, xCol, yCol, plotTitleStr, plotNameStr, dataSourceStr,
                         saveDir, bandList, threshStruct, instrument, plotWidth=280, nBins=50):
    if xCol == "detector":
        visitDataAll = visitDataAll.astype({xCol: int})
    if yCol == "detector":
        visitDataAll = visitDataAll.astype({yCol: int})

    visitList = list(set(visitDataAll["visit"]))
    if len(str(visitList[0])) == 13:
        obsDateFirst = str(min(visitList))[:8]
        obsDateLast = str(max(visitList))[:8]
    else:
        obsDateFirst = None
        obsDateLast = None

    color = "teal"
    lineColors = ["black", "white"]
    bandColorDict = setDesBandColorDict()

    axisScaleDict = setAxisScaleDict()
    xAxScale = axisScaleDict[xCol]
    yAxScale = axisScaleDict[yCol]
    isAstrom = "astromOffset" in xCol and "astromOffset" in yCol
    gridList = []
    # Make some reasonable limits common to all plots (exceptions in loop)
    xDataMin = np.nanmin(visitDataAll[xCol])
    xDataMax = np.nanmax(visitDataAll[xCol])
    xDataRange = np.abs(xDataMin - xDataMax)
    yDataMin = np.nanmin(visitDataAll[yCol])
    yDataMax = np.nanmax(visitDataAll[yCol])
    yDataRange = np.abs(yDataMin - yDataMax)
    rangeScale = 0.2
    if xAxScale == "log":
        xDataMin *= (1 - rangeScale)
        xDataMax *= (1 + rangeScale)
    else:
        xDataMin -= rangeScale*xDataRange
        xDataMax += rangeScale*xDataRange
    if yAxScale == "log":
        yDataMin *= (1 - rangeScale)
        yDataMax *= (1 + rangeScale)
    else:
        yDataMin -= rangeScale*yDataRange
        yDataMax += rangeScale*yDataRange

    rayleighConst = np.sqrt(np.pi/2)/np.sqrt((2 - np.pi/2))
    for ib, band in enumerate(bandList):
        visitBandAll = visitDataAll[visitDataAll["band"] == band]
        visitBand = visitBandAll[np.isfinite(visitBandAll[xCol])].copy()

        source = getSource(visitBand, [xCol], [yCol])

        nPoints = 1000 if isAstrom else 10
        # Let these guys have per-band limits instead of sharing with all bands
        if xCol in ["zeroPoint", "skyBg", "expTimeScaledSkyBg"]:
            xDataMin = np.nanmin(visitBand[xCol])
            xDataMax = np.nanmax(visitBand[xCol])
            xDataRange = np.abs(xDataMin - xDataMax)
            xDataMin -= rangeScale*xDataRange
            xDataMax += rangeScale*xDataRange
        if yCol in ["zeroPoint", "skyBg", "expTimeScaledSkyBg="]:
            yDataMin = np.nanmin(visitBand[yCol])
            yDataMax = np.nanmax(visitBand[yCol])
            yDataRange = np.abs(yDataMin - yDataMax)
            yDataMin -= rangeScale*yDataRange
            yDataMax += rangeScale*yDataRange
        # Make histogram for right plot
        histCol = yCol if xCol in ["visit", "detector"] else xCol
        visitBand = visitBand[np.isfinite(visitBand[histCol])]
        xHistMin = xDataMin if histCol == xCol else yDataMin
        xHistMax = xDataMax if histCol == xCol else yDataMax
        histScale = axisScaleDict[histCol]
        if histScale == "log":
            xHistMin, xHistMax = np.log10(xHistMin), np.log10(xHistMax)
        histData = np.log10(visitBand[histCol]) if histScale == "log" else visitBand[histCol]
        hist, edges = np.histogram(histData, density=False, bins=nBins)

        unitsDict = setUnitsDict()
        xUnitFormat = "%.2f"
        if unitsDict[xCol] in ["number", "degrees", "counts"]:
            xUnitFormat = "%.1e" if xCol == "visit" else "%d"
        if xCol in ["psfStarScaledDeltaSizeScatter", "medianE"] or "astrom" in xCol:
            xUnitFormat = "%.3f"
        yUnitFormat = "%.2f"
        if unitsDict[yCol] in ["number", "degrees", "counts"]:
            yUnitFormat = "%.1e" if yCol == "visit" else "%d"
        if yCol in ["psfStarScaledDeltaSizeScatter", "medianE"] or "astrom" in yCol:
            yUnitFormat = "%.3f"

        xList = np.linspace(xDataMin, xDataMax, 500)
        yMean = [np.nanmean(source.data[yCol])]*len(xList)
        yStd = [np.nanstd(source.data[yCol])]*len(xList)
        yStdPlus = [x + y for x, y in zip(yMean, yStd)]
        yStdMinus = [x - y for x, y in zip(yMean, yStd)]

        if unitsDict[yCol] in ["pixels", "degrees", "counts"] or "effTime" in yCol:
            meanLabel = f"mean = {yUnitFormat}" % yMean[0]
            stdLabel = f"std = {yUnitFormat}" % yStd[0]
        else:
            unitFact = 1000.0  # convert to units of "milli"
            yUnitStr = "as" if unitsDict[yCol] == "arcsec" else unitsDict[yCol]
            if yUnitStr == "":
                yUnitStr = "illi"
            meanLabel = f"mean = {yUnitFormat}" % (yMean[0]*unitFact)
            meanLabel += f" (m{yUnitStr})"
            stdLabel = f"std = {yUnitFormat}" % (yStd[0]*unitFact)
            stdLabel += f" (m{yUnitStr})"

        if isAstrom:
            yExpected = (1/rayleighConst)*xList
            xExpectedLog = np.log10(xList)
            binWidth = edges[1] - edges[0]
            sigma = np.nanmean(source.data[xCol])/np.sqrt(np.pi/2)
            pdfRayleigh = (xList/sigma**2)*np.exp(-xList**2/(2*sigma**2))
            pdfRayleigh *= len(source.data[xCol])*binWidth/len(edges)
            measuredSlope = np.nanmean(source.data[xCol]/source.data[yCol])
            meanLabel = "meas = {:.2f}".format(measuredSlope)
            yMean = 1/measuredSlope*xList
            yHistPlot = [0, 1.05*max(hist)]
            xMean = [np.nanmean(source.data[xCol])]*len(yHistPlot)
            xMeanLog = np.log10(xMean)
        # create left scatter plot
        color = bandColorDict[band]
        # color = palettes.Colorblind[7][ib]
        nonSelectColor = palettes.Pastel1[7][ib]
        left = figure(width=plotWidth, height=plotWidth,
                      x_axis_type=xAxScale, y_axis_type=yAxScale,
                      title="band: {}    N = {}".format(band, len(source.data[xCol])))
        left.circle(xCol, yCol, source=source, fill_alpha=0.2, color=color, size=6,
                    hover_color=HOVER_COLOR,
                    selection_fill_color=color, # selection_line_color=color,
                    nonselection_fill_color=nonSelectColor,
                    nonselection_line_color=nonSelectColor,
                    name="toHover")
        left.x_range = Range1d(xDataMin, xDataMax)
        left.y_range = Range1d(yDataMin, yDataMax)
        if isAstrom:
            left.line(xList, yExpected, line_color=lineColors[0], line_width=1,
                      line_dash="1 5",
                      legend_label="exp = {:.2f}".format(rayleighConst))
        left.line(xList, yMean, alpha=1, line_color="black", line_width=1,
                  line_dash=(6, 6), legend_label=meanLabel)
        if not isAstrom:
            left.line(xList, yStdPlus, alpha=1, line_color="black", line_width=1,
                     line_dash=(1, 5), legend_label=stdLabel)
            left.line(xList, yStdMinus, alpha=1, line_color="black", line_width=1,
                     line_dash=(1, 5), legend_label=stdLabel)
        left, xNumRej, yNumRej, totNumRej = addBokehThresholds(
            xCol, yCol, left, threshStruct, xData=visitBand[xCol], yData=visitBand[yCol]
        )

        TOOLTIPS_LIST = [("visit", "@visit"), ("detector", "@detector"), ("band", "@band")]
        hoverLeft = HoverTool(name="toHover", tooltips=TOOLTIPS_LIST)
        left.add_tools(hoverLeft)
        # left.tools.append(hoverLeft)
        left.legend.location = "top_left"
        left.legend.padding = 0
        left.legend.spacing = 0
        if unitsDict[xCol] in [""]:
            left.xaxis.axis_label = "{} [{}]".format(xCol, band)
        else:
            left.xaxis.axis_label = "{} ({}) [{}]".format(xCol, unitsDict[xCol], band)
        if unitsDict[yCol] in [""]:
            left.yaxis.axis_label = "{} [{}]".format(yCol, band)
        else:
            left.yaxis.axis_label = "{} ({}) [{}]".format(yCol, unitsDict[yCol], band)
        left.xaxis[0].formatter = PrintfTickFormatter(format=xUnitFormat)
        left.yaxis[0].formatter = PrintfTickFormatter(format=yUnitFormat)
        left.xaxis.major_label_orientation = np.pi/4

        # create right histogram
        if ib%2 == 0 :
            rTitle = plotTitleStr
        else:
            if obsDateFirst is not None:
                if obsDateFirst == obsDateLast:
                    rTitle = "dayObs: {}".format(obsDateFirst)
                else:
                    rTitle = "dayObs: {}_{}".format(obsDateFirst, obsDateLast)
            else:
                rTitile = ""
        right = figure(width=plotWidth, height=plotWidth,
                       title="{}".format(rTitle))
        right.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
                   fill_color=color, line_color=color,
                   hover_color=HOVER_COLOR, hover_alpha=0.8,
                   nonselection_fill_color=nonSelectColor,
                   nonselection_line_color=nonSelectColor)
        right.x_range = Range1d(xHistMin, xHistMax)
        if isAstrom:
            if 0:  # This is not working right...
                right.line(xExpectedLog, pdfRayleigh, line_color="black", line_width=1,
                           line_dash="1 5")
            right.line(xMean, yHistPlot, line_color="black", line_width=1,
                       line_dash="6 6", legend_label="mean")
        right, xNumRej, yNumRej, totNumRej = addBokehThresholds(
            xCol, yCol, right, threshStruct, xData=visitBand[xCol], yData=None
        )
        if xNumRej is not None:
            xStart = right.x_range.start
            yStart = 0.0
            xDelta = right.x_range.end - xStart
            yDelta = max(hist)
            if not np.isfinite(xDelta):
                xDelta = 1.05*(np.nanmax(visitBand[xCol]) - np.nanmin(visitBand[xCol]))
                xStart = np.nanmin(visitBand[xCol])
            rejColor = "royalblue"
            nRejKwargs = dict(x_units="data", y_units="data",
                              text_font_size="12px", text_color="white",
                              border_line_color=rejColor, border_line_alpha=1.0,
                              background_fill_color=rejColor, background_fill_alpha=1.0,
                              border_line_width=6,
                              )
            xNumRejLabel = Label(
                x=xStart + 0.97*xDelta, y=yStart + 0.9*yDelta,
                text="numRej = {}".format(xNumRej), text_align="right", **nRejKwargs
            )
            right.add_layout(xNumRejLabel)

        unitFormat = "%.2f"
        if unitsDict[histCol] in ["number", "degrees", "counts"]:
            unitFormat = "%.1e" if histCol == "visit" else "%d"
        right.xaxis[0].formatter = PrintfTickFormatter(format=unitFormat)
        right.yaxis[0].formatter = PrintfTickFormatter(format="%d")
        rAxisLabel = ("log10({}) ({}) [{}]".format(histCol, unitsDict[histCol], band) if histScale == "log"
                      else "{} ({}) [{}]".format(histCol, unitsDict[histCol], band))
        right.xaxis.axis_label = rAxisLabel
        right.yaxis.axis_label = "Number"
        right.xaxis.major_label_orientation = np.pi/4
        if ib == 0:
            gridList.append([left, right])
        elif np.mod(ib, 2) > 0:
            gridList[-1].extend([left, right])
        else:
            gridList.append([left, right])
    p = gridplot(gridList)

    subSaveDir = "histogramPlots/"
    if not os.path.exists(saveDir + subSaveDir):
        print("Creating new directory: {}".format(saveDir + subSaveDir))
        os.mkdir(saveDir + subSaveDir)

    plotFilename = "{}{}{}_{}_{}_vs_{}_{}.html".format(
        saveDir, subSaveDir, instrument, dataSourceStr, yCol, xCol, plotNameStr)
    output_file(plotFilename,
                title="{} {} {}: {} vs. {}".format(instrument, plotNameStr, dataSourceStr, xCol, yCol),
                mode="inline")
    print("Saving {}".format(plotFilename))
    save(p)
    return None


def makeCorrPlots(visitDataAll, xColList, yColList, plotTitleStr, plotNameStr, dataSourceStr,
                  saveDir, instrument, butler=None, xColHistList=None, yColHistList=None,
                  setXRange=False, setYRange=False, collection=None
):
    nDataId = len(visitDataAll)
    # fileExtensionList = ["png"] if nDataId > N_DATA_HTML_MAX else ["html", "png"]
    # fileExtensionList = ["png"] if nDataId > N_DATA_HTML_MAX else ["html"]
    fileExtensionList = ["html", "png"]
    threshStruct = introspectConfigs(instrument, butler=None)
    print(threshStruct)
    # Make a sorted list of the bands that actually exist
    bandList = list(set(visitDataAll["band"]))
    bandList = [band for band in SORTED_FULL_BAND_LIST if band in bandList]
    print("Sorted list of existing bands: {}".format(bandList))
    for fileExtension in fileExtensionList:
        makeParameterCompBokehPlots(
            visitDataAll, fileExtension, xColList, yColList,
            plotTitleStr, plotNameStr, dataSourceStr,
            saveDir, bandList, threshStruct, instrument,
            collection=collection, setXRange=setXRange, setYRange=setYRange
        )
    if xColHistList is not None and yColHistList is not None:
        for xCol, yCol in zip(xColHistList, yColHistList):
            makeScatAndHistPlots(
                visitDataAll, xCol, yCol,
                plotTitleStr, plotNameStr, dataSourceStr,
                saveDir, bandList, threshStruct, instrument
            )
    return None


def makeHistPlots(visitDataAll, xColList, yColList, plotTitleStr, plotNameStr, dataSourceStr,
                  saveDir, instrument, butler=None, xColHistList=None, yColHistList=None):
    threshStruct = introspectConfigs(instrument, butler=None)
    # Make a sorted list of the bands that actually exist
    bandList = list(set(visitDataAll["band"]))
    bandList = [band for band in SORTED_FULL_BAND_LIST if band in bandList]
    if xColHistList is not None and yColHistList is not None:
        for xCol, yCol in zip(xColHistList, yColHistList):
            makeScatAndHistPlots(
                visitDataAll, xCol, yCol,
                plotTitleStr, plotNameStr, dataSourceStr,
                saveDir, bandList, threshStruct, instrument
            )
    return None
