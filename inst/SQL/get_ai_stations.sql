-- Query previously successful stations in the Aleutian Isl. Bottom Trawl Survey
-- Authors: Ned Laman, Zack Oyafuso

(
-- First, query trawlable stations with well-performing hauls from RACEBASE.HAUL
-- that occurred 1994-current. Remove hauls that have a NULL value for
-- stationid and stratum. The 1994 cutoff seems to be related to a change from
-- LORAN in 1991 to GPS in 1994.
 SELECT
 DISTINCT STATIONID, STRATUM
 FROM RACEBASE.HAUL
 WHERE REGION = 'AI'
 AND CRUISE >= 199401
 AND PERFORMANCE >= 0
 AND STATIONID IS NOT NULL
 AND STRATUM IS NOT NULL
 AND (STATIONID, STRATUM) IN (SELECT DISTINCT STATIONID, TO_NUMBER(STRATUM) STRATUM FROM AI.AIGRID_GIS)
 GROUP BY STATIONID, STRATUM

UNION

 -- Unite with trawlable stations with well-performing hauls from 1991 using
 -- the Ocean Hope that have yet to be replicated 1994-on.
 SELECT
 STATIONID,
 STRATUM
 FROM AI.AIGRID_GIS
 WHERE TRAWLABLE = 'Y'
 AND (STRATUM, STATIONID) IN
  (SELECT
   DISTINCT STRATUM,
   STATIONID
   FROM RACEBASE.HAUL
   WHERE REGION = 'AI'
   AND PERFORMANCE >= 0
   AND CRUISE = 199101
   AND HAUL_TYPE = 3
   AND STATIONID IS NOT NULL
   AND STRATUM IS NOT NULL
   AND VESSEL = 85

   MINUS

   SELECT
   DISTINCT STRATUM,
   STATIONID
   FROM RACEBASE.HAUL
   WHERE PERFORMANCE >= 0
   AND HAUL_TYPE = 3
   AND STRATUM IS NOT NULL
   AND STATIONID IS NOT NULL
   AND CRUISE > 199400
   AND REGION = 'AI'
  )
)

-- Remove untrawlable stations from this concatenated list of stations.
MINUS

SELECT
DISTINCT STATIONID,
STRATUM
FROM AI.AIGRID_GIS
WHERE TRAWLABLE = 'N'

ORDER BY STRATUM, STATIONID
