#!/usr/bin/env python
import argparse
import asyncio
import bokehPlotUtils as bkPlotUtils
from datetime import datetime, timezone, timedelta

import os
import sys

from lsst_efd_client import EfdClient

sys.path.append(os.path.abspath("~/Python"))

# import lsst.log as lsstLog
# logger = lsstLog.Log.getLogger("makeBokehVisitParCorrPlot.py")


async def readFromSasquatch(efdClientName, dataBaseName, dataSourceStr, datetimeStart, datetimeEnd):
    dataType = f"{dataBaseName}.{dataSourceStr}"
    efd_client = EfdClient(efdClientName, db_name=dataBaseName)
    query = f"SELECT * FROM \"{dataType}\" WHERE time > '{datetimeStart}' AND " \
            f"time < '{datetimeEnd}' AND dataset_tag='LSSTComCamSim/rapid_analysis'"
    print("Reading data from Sasquatch with query:")
    print(query)
    result = await efd_client.influx_client.query(query)
    return result


def readFromConsDb(dayObsStart, dayObsEnd):
    from lsst.summit.utils import ConsDbClient
    with open("/sdf/home/l/laurenma/.consDB_token", "r") as f:
        token = f.read()
    consDbClient = ConsDbClient(f"https://user:{token}@usdf-rsp.slac.stanford.edu/consdb")

    print(consDbClient.schema())
    dayObsStartInt = int(dayObsStart.replace("-", ""))
    dayObsEndInt = int(dayObsEnd.replace("-", ""))
    instrument = "lsstcomcamsim"
    whereStr = f'''
        WHERE obs_start_mjd IS NOT NULL
        AND s_ra IS NOT NULL
        AND s_dec IS NOT NULL
        AND sky_rotation IS NOT NULL
        AND ((band IS NOT NULL) OR (physical_filter IS NOT NULL))
        AND day_obs >= {dayObsStartInt}
    '''

    exposure_query = f'''SELECT * FROM cdb_{instrument}.exposure {whereStr}'''
    visit1_query = f'''SELECT * FROM cdb_{instrument}.visit1 {whereStr}'''
    ccdexposure_quicklook_query = f'''SELECT * FROM cdb_{instrument}.ccdexposure_quicklook'''
    ccdvisit1_quicklook_query = f'''SELECT * FROM cdb_{instrument}.ccdvisit1_quicklook'''


def main(dayObsStart=None, datetimeStart=None, datetimeEnd=None, instrument=None,
         repo=None, collections=None, skymapName=None,
         # readFromButler=False,
         saveDir=None, doVisitSkyMapPlots=False):

    currentUtc = datetime.now(timezone.utc)
    currentYear = currentUtc.year
    if datetimeStart is None:
        if dayObsStart is None:  # Set to current UTC day if not set explicitly
            dayObsStart = str(currentUtc.date())
        datetimeStart = datetime.fromisoformat(str(dayObsStart))
    else:
        dayObsStart = datetimeStart[:10]

    if datetimeEnd is None or datetimeEnd == "None":
        datetimeEnd = datetimeStart + timedelta(days=1)
    dayObsEnd = str(datetimeEnd)[:10]

    print("datetimeStart = {} datetimeEnd = {}".format(datetimeStart, datetimeEnd))
    print("  dayObsStart = {}            dayObsEnd = {}".format(dayObsStart, dayObsEnd))

    efdClientName = "usdfdev_efd"
    dataBaseName = "lsst.dm"
    dataSourceStr = "calexpSummaryMetrics"

    # if not readFromButler and (datetimeStart is None or datetimeStart == "None"):
    if datetimeStart is None or datetimeStart == "None" or datetimeEnd is None or datetimeEnd == "None":
        raise RuntimeError("Must specify either --dayObsStart or --datetimeStart AND --dataStartEnd")
    # if readFromButler and (butler is None or collections is None):
    #     raise RuntimeError("Must provide --butler and --collections if --readFromButler is True")

    visitDataAll = asyncio.run(
        readFromSasquatch(efdClientName, dataBaseName, dataSourceStr, datetimeStart, datetimeEnd)
    )
    if len(visitDataAll) == 0:
        raise RuntimeError("No data found in Sasquatch with search constraints.")

    # Don't assume data are sorted by visitId
    visitDataAll.sort_values(["visit", "detector"], ascending=[True, True], inplace=True)

    print("dayObsStart = {}  dayObsEnd = {}".format(dayObsStart, dayObsEnd))
    dayObsStartStrip = dayObsStart.replace("-", "")
    dayObsEndStrip = dayObsEnd.replace("-", "")
    # LSSTComCamSim data begin with 7024 (i.e. not 2024), so we need to account
    # for the this...
    visitList = list(set(visitDataAll["visit"]))
    minVisitYear = int(str(min(visitList))[:4])
    yearOffset = int((minVisitYear - currentYear)/1000)
    visitDayObsStartInt = int(dayObsStartStrip) + yearOffset*10000000
    visitDayObsEndInt = int(dayObsEndStrip) + yearOffset*10000000
    visitSelectMask = []
    for visit in visitDataAll["visit"].values:
        visitDayObsInt = int(visit[:8])
        if visitDayObsInt >= visitDayObsStartInt and visitDayObsInt <= visitDayObsEndInt:
            visitSelectMask.append(True)
        else:
            visitSelectMask.append(False)
    visitDataAll = visitDataAll[visitSelectMask].copy()
    if len(visitDataAll) == 0:
        raise RuntimeError("No data found with search constraints.")
    visitDataAll = bkPlotUtils.addMetricColumns(visitDataAll)
    visitDataAll = bkPlotUtils.addSeqNumColumn(visitDataAll)
    if instrument is not None:
        print("Selecting data from instrument: {}".format(instrument))
        visitDataAll = visitDataAll[visitDataAll["instrument"].values == instrument]
    else:
        print("WARNING: no instrument has been specified.  Use --instrument INSTRUMENTNAME to set one.")
    nDataId = len(visitDataAll)
    if nDataId == 0:
        raise RuntimeError("No data found with search and instrument constraints.")

    visitList = list(set(visitDataAll["visit"]))
    print("minVisit = {}  maxVisit = {}".format(min(visitList), max(visitList)))
    if not doVisitSkyMapPlots:
        plotTitleStr = "Per-Detector effTime Metrics"
        plotNameStr = "effTimeMetrics" + "_" + dayObsStart + "_" + dayObsEnd
        xColList = ["seqNum", "seqNum", "seqNum",
                    "seqNum", "seqNum", "seqNum",
                    "seqNum", "seqNum"]
        yColList = ["psfFwhm", "expTimeScaledZp", "skyBg",
                    "effTimePsfSigmaScale", "effTimeZeroPointScale", "effTimeSkyBgScale",
                    "effTime", "effTimeScale"]
        bkPlotUtils.makeCorrPlots(
            visitDataAll, xColList, yColList, plotTitleStr, plotNameStr, dataSourceStr, saveDir,
            instrument=instrument, butler=None, setXRange=True, setYRange=False
        )

        plotTitleStr = "Per-Detector Metric Correlations"
        plotNameStr = "metricCorrelations" + "_" + dayObsStart + "_" + dayObsEnd
        xColList = ["psfApFluxDelta", "psfApFluxDelta", "psfApFluxDelta",
                    "psfApCorrSigmaScaledDelta", "psfApCorrSigmaScaledDelta", "psfApCorrSigmaScaledDelta",
                    "psfApFluxDelta"]
        yColList = ["medianE", "psfStarScaledDeltaSizeScatter", "psfFwhm",
                    "medianE", "psfStarScaledDeltaSizeScatter", "psfFwhm",
                    "psfApCorrSigmaScaledDelta"]
        # xColList = ["astromOffsetMean", "zenithDistance", "skyBg",
        #             "zenithDistance", "detector", "seqNum",
        #             "psfFwhm", "psfFwhm", "psfTraceRadiusDelta"]
        # yColList = ["astromOffsetStd", "astromOffsetMean", "astromOffsetMean",
        #             "psfFwhm", "expTimeScaledZp", "expTimeScaledZp",
        #             "medianE", "psfStarScaledDeltaSizeScatter", "maxDistToNearestPsf"]
        bkPlotUtils.makeCorrPlots(
            visitDataAll, xColList, yColList, plotTitleStr, plotNameStr, dataSourceStr, saveDir,
            instrument=instrument, butler=None, setXRange=True, setYRange=True
        )

        plotTitleStr = "Per-Detector Metric-Parameter Correlations"
        plotNameStr = "parameterMetricCorrelations" + "_" + dayObsStart + "_" + dayObsEnd
        yColList = ["psfApFluxDelta", "psfApCorrSigmaScaledDelta", "psfStarScaledDeltaSizeScatter",
                    "medianE", "psfTraceRadiusDelta", "astromOffsetMean",
                    "psfApFluxDelta", "psfApCorrSigmaScaledDelta", "psfStarScaledDeltaSizeScatter",
                    "medianE", "psfTraceRadiusDelta", "astromOffsetMean",
                    "psfApFluxDelta", "psfApCorrSigmaScaledDelta", "psfStarScaledDeltaSizeScatter",
                    "medianE", "psfTraceRadiusDelta", "astromOffsetMean",
        ]
        xColList = ["psfFwhm", "psfFwhm", "psfFwhm",
                    "psfFwhm", "psfFwhm", "psfFwhm",
                    "zenithDistance", "zenithDistance", "zenithDistance",
                    "zenithDistance", "zenithDistance", "zenithDistance",
                    "skyBg", "skyBg", "skyBg",
                    "skyBg", "skyBg", "skyBg",
        ]
        bkPlotUtils.makeCorrPlots(
            visitDataAll, xColList, yColList, plotTitleStr, plotNameStr, dataSourceStr, saveDir,
            instrument=instrument, butler=None, setXRange=True, setYRange=False
        )

        plotTitleStr = "Per-Detector Parameter Correlations"
        plotNameStr = "parameterCorrelations" + "_" + dayObsStart + "_" + dayObsEnd
        xColList = ["astromOffsetMean", "zenithDistance", "skyBg",
                    "zenithDistance", "detector", "seqNum",
                    # "psfFwhm", "psfFwhm", "psfTraceRadiusDelta"
        ]
        yColList = ["astromOffsetStd", "astromOffsetMean", "astromOffsetMean",
                    "psfFwhm", "expTimeScaledZp", "expTimeScaledZp",
                    # "medianE", "psfStarScaledDeltaSizeScatter", "maxDistToNearestPsf"
        ]
        bkPlotUtils.makeCorrPlots(
            visitDataAll, xColList, yColList, plotTitleStr, plotNameStr, dataSourceStr, saveDir,
            instrument=instrument, butler=None, setXRange=False, setYRange=False
        )


        # List of parameters to make histogram plots for:
        plotTitleStr = "Per-Band Per-Detector Metrics"
        plotNameStr = dayObsStart + "_" + dayObsEnd
        xColHistList = ["zenithDistance", "astromOffsetMean", "skyBg",
                        "zenithDistance", "detector", "detector",
                        "psfFwhm", "psfFwhm"]
        yColHistList = ["astromOffsetMean", "astromOffsetStd", "astromOffsetMean",
                        "psfFwhm", "effTimeZeroPointScale", "effTimePsfSigmaScale",
                        "medianE", "psfStarScaledDeltaSizeScatter"]

        bkPlotUtils.makeHistPlots(
            visitDataAll, xColList, yColList, plotTitleStr, plotNameStr, dataSourceStr, saveDir,
            instrument=instrument, butler=None, xColHistList=xColHistList, yColHistList=yColHistList
        )

    if doVisitSkyMapPlots:
        # Make showVisitSkyMap plots for current list of available calexps.
        if collections is None:
            from lsst.daf.butler import Butler
            from lsst.daf.butler.registry import MissingCollectionError

            butler = Butler(repo)
            collections = list(set(visitDataAll["run"]))
            collectionsStr = ""
            for collection in collections:
                try:
                    butler.registry.getCollectionType(collection)
                    collectionsStr = collectionsStr + " " + str(collection)
                except MissingCollectionError:
                    print("No collection with name {} found in repo: {}.  "
                          "Not adding it to the collections list.".format(collection, repo))
                    missingCollectionData = visitDataAll[visitDataAll["run"] == collection].copy(deep=True)
                    omitVisitList = list(set(missingCollectionData["visit"]))
                    print("omitVisitList = ", omitVisitList)
                    for visit in omitVisitList:
                        if visit in visitList:
                            visitList.remove(visit)
        else:
            collectionsStr = collections
        # Might only have RUN collections without camera collection CHAINED in.
        # We can usually count on there being an INSTRUMENT/defaults collection
        # being defined, so add it to the collection list.
        defaultsCollection = instrument + "/defaults"
        try:
            butler.registry.getCollectionType(defaultsCollection)
            collectionsStr = collectionsStr + " " + defaultsCollection
        except MissingCollectionError:
            print("No defaults collection with name {} found in repo: {}.  "
                  "Not adding it to the collections list (so showVisitSummay may fail to find a camera.".
                  format(defaultsCollection, repo))

        try:
            skyButler = Butler(repo, collections="skymaps", instrument=instrument, skymap=skymapName)
            skyMap = skyButler.get("skyMap", instrument=instrument, skymap=skymapName)
        except MissingCollectionError:
             skyMap = None
             doPatchOutline = True
             print("Unable to locate skymap. Plotting patch outlines...")
        if skyMap is not None:
            patchSizePix = max(
                skyMap[0].getPatchInfo(0).outer_bbox.width, skyMap[0].getPatchInfo(0).outer_bbox.height
            )
            patchSizeDeg = skyMap.config.pixelScale*patchSizePix/3600.0
            doPatchOutline = True if 15*patchSizeDeg < 1.8 else False

        SKYMAP_DIR = os.environ["SKYMAP_DIR"]
        # Make single plot for all visits
        visitsStr = ""
        for visit in visitList:
            visitsStr = visitsStr + " " + str(visit)

        print("collectionsStr = {}".format(collectionsStr))
        plotNameStr = f"showVisit_{instrument}_{dayObsStart}_{dayObsEnd}"
        saveFile = f" --saveFile {saveDir}{plotNameStr}.png"
        args = f" {repo} --collections {collectionsStr} --skymapName {skymapName}"
        args += f" --visits{visitsStr}"  # --doUnscaledLimitRatio"
        args += saveFile
        cmd = f"{SKYMAP_DIR}/doc/_static/skymap/showVisitSkyMap.py"
        print("Running {} {}".format(cmd, args))
        os.system(cmd + args)

        # Subdivide plots in to RA/Dec groupings
        maxDeltaDeg = 6.0
        ras = [int(ra) for ra in visitDataAll["ra"].values]
        decs = [int(dec) for dec in visitDataAll["dec"].values]
        raDecList = []
        for ra, dec in zip(ras, decs):
            if [ra, dec] not in raDecList:
                raDecList.append([ra, dec])
        raDecList.sort()
        raDecListTrimmed = [raDecList[0]]
        index = 0
        for raDec in raDecList[1:]:
            if raDec[0] < 360 - maxDeltaDeg/2 and raDec[0] > maxDeltaDeg/2 :
                if (abs(raDec[0] - raDecListTrimmed[index][0]) > maxDeltaDeg
                    or abs(raDec[1] - raDecListTrimmed[index][1]) > maxDeltaDeg):
                    raDecListTrimmed.append(raDec)
                    index += 1
            else:
                if abs(raDec[1] - raDecListTrimmed[index][1]) > maxDeltaDeg:
                    raDecListTrimmed.append(raDec)
                    index += 1

        print("Found {} RA/Dec subgroups: {}".format(len(raDecListTrimmed), raDecListTrimmed))
        raDecDict = {}
        for raDec in raDecListTrimmed:
            raDecKey = str(raDec)
            raRef = raDec[0]
            decRef = raDec[1]
            raDecDict[raDecKey] = []
            for index, row in visitDataAll.iterrows():
                visit = row["visit"]
                ra = row["ra"]
                dec = row["dec"]
                if visit not in raDecDict[raDecKey]:
                    if raRef < 360 - maxDeltaDeg/2 and raRef > maxDeltaDeg/2 :
                        if abs(ra - raRef) < maxDeltaDeg and abs(dec - decRef) < maxDeltaDeg:
                            raDecDict[raDecKey].append(visit)
                    else:
                        if abs(dec - decRef) < maxDeltaDeg:
                            raDecDict[raDecKey].append(visit)

        for i, (raDecKey, subVisitList) in enumerate(raDecDict.items()):
            visitsStr = ""
            for visit in subVisitList:
                visitsStr = visitsStr + " " + str(visit)

            plotNameStr = f"showVisit_{instrument}_{dayObsStart}_{dayObsEnd}"
            plotNameStr += "_" + str(i)
            saveFile = f" --saveFile {saveDir}{plotNameStr}.png"
            args = f" {repo} --collections {collectionsStr} --skymapName {skymapName}"
            args += f" --visits{visitsStr}"
            if doPatchOutline:
                args += " --showPatch"
            args += saveFile
            # cmd = "/sdf/home/l/laurenma/Python/showVisitSkyMap.py"
            cmd = f"{SKYMAP_DIR}/doc/_static/skymap/showVisitSkyMap.py"
            print("Running {} {}".format(cmd, args))
            os.system(cmd + args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dayObsStart", type=str, default=None,
                        help=("UTC date for data search.  Format: \"YYYY-MM-DD\"\n"
                              "Note this is the date of data upload to Sasquatch (which could have a "
                              " delay from the observation date)."))
    parser.add_argument("--datetimeStart", type=str, default=None,
                        help=("Start UTC datetime for data search.  Format: \"YYYY-MM-DD HH:MM:SS\"\n"
                              "Note this is the date of data upload to Sasquatch (which could have a "
                              " delay from the observation date)."))
    parser.add_argument("--datetimeEnd", type=str, default="2024-06-30 23:59:59",
                        help=("End UTC datetime for data search.  Format: \"YYYY-MM-DD HH:MM:SS\"\n"
                              "Note this is the date of data upload to Sasquatch (which could have a "
                              " delay from the observation date)."))
    parser.add_argument("--instrument", type=str, default="LSSTComCamSim",
                        help="Name of the instrument for the data search.")
    parser.add_argument("--repo", type=str, default="embargo_or4",  # "/repo/embargo",
                        help="URI or path to an existing data repository root or configuration file")
    parser.add_argument("--collections", type=str, nargs="+", default=None,
                        help="Blank-space separated list of collection names for butler instantiation",
                        metavar=("COLLECTION1", "COLLECTION2"), required=False)
    parser.add_argument("--skymapName", default="ops_rehersal_prep_2k_v1",
                        help="Name of the skymap for the collection")
    # parser.add_argument("--readFromButler", action="store_true", default=False,
    #                     help="Read from a butler repo instead of Sasquatch if True.
    #                     TODO: add support for this if deemed worth doing.")
    parser.add_argument("--saveDir", type=str,
                        default="/sdf/home/l/laurenma/public_html/OpsRehearsal4/monitor/",
                        help=("Full path to directory for saving plots"))
    parser.add_argument("--doVisitSkyMapPlots", action="store_true", default=False,
                        help="Make showVisitSkyMap plots if True.  Otherwise Bokeh plots are made.")
    args = parser.parse_args()
    main(dayObsStart=args.dayObsStart, datetimeStart=args.datetimeStart, datetimeEnd=args.datetimeEnd,
         instrument=args.instrument, repo=args.repo, collections=args.collections,
         skymapName=args.skymapName,  # readFromButler=args.readFromButler,
         saveDir=args.saveDir, doVisitSkyMapPlots=args.doVisitSkyMapPlots)
