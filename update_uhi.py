import ee, json, os, datetime as dt

SERVICE_ACCOUNT = "uhi-updater@nsti-513504.iam.gserviceaccount.com"
PROJECT = "nsti-513504"

# West, South, East, North
AREA_A = [55.2004, 25.1114, 55.2065, 25.1149]   # Target (Al Barsha 1)
AREA_R = [55.1997, 25.1001, 55.2058, 25.1038]   # Reference

START = "2025-06-01"
HAZE_MIN_C = 40               # drop likely cloud/haze scenes below this (hot months only)
HAZE_MONTHS = (5, 6, 7, 8, 9)

creds = ee.ServiceAccountCredentials(SERVICE_ACCOUNT, key_data=os.environ["GEE_KEY"])
ee.Initialize(creds, project=PROJECT)

area_a = ee.Geometry.Rectangle(AREA_A)
area_r = ee.Geometry.Rectangle(AREA_R)
end = (dt.date.today() + dt.timedelta(days=1)).isoformat()

def lst_image(img):
    qa = img.select("QA_PIXEL")
    clear = (qa.bitwiseAnd(1).eq(0)       # no fill
             .And(qa.bitwiseAnd(8).eq(0))   # no cloud
             .And(qa.bitwiseAnd(16).eq(0))) # no cloud shadow
    st = img.select("ST_B10")
    lst = st.multiply(0.00341802).add(149.0).subtract(273.15).rename("LST_C")
    return lst.updateMask(clear.And(st.gt(0))).copyProperties(img, ["system:time_start"])

l9 = (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
      .filterBounds(area_a).filterDate(START, end).map(lst_image))

def per_scene(img):
    a = img.reduceRegion(ee.Reducer.mean(), area_a, 30, maxPixels=1e9).get("LST_C")
    r = img.reduceRegion(ee.Reducer.mean(), area_r, 30, maxPixels=1e9).get("LST_C")
    return ee.Feature(None, {"date": img.date().format("YYYY-MM-dd"), "a": a, "r": r})

fc = l9.map(per_scene).filter(ee.Filter.notNull(["a", "r"]))
raw = sorted([f["properties"] for f in fc.getInfo()["features"]], key=lambda p: p["date"])

scenes = []
for p in raw:
    month = int(p["date"][5:7])
    if month in HAZE_MONTHS and p["r"] < HAZE_MIN_C:
        continue
    scenes.append({
        "date": p["date"],
        "area_a_c": round(p["a"], 2),
        "area_r_c": round(p["r"], 2),
        "difference_c": round(p["a"] - p["r"], 2),
    })

# NDVI from Sentinel-2 (median composite)
s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(area_a).filterDate(START, end)
      .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
      .map(lambda i: i.normalizedDifference(["B8", "B4"]).rename("NDVI")))
ndvi = s2.median()

def mean_ndvi(geom):
    v = ndvi.reduceRegion(ee.Reducer.mean(), geom, 10, maxPixels=1e9).get("NDVI").getInfo()
    return round(v, 3) if v is not None else None

def avg(key):
    return round(sum(s[key] for s in scenes) / len(scenes), 2) if scenes else None

out = {
    "updated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    "period_start": START,
    "latest_scene_date": scenes[-1]["date"] if scenes else None,
    "scenes_available": len(raw),
    "scenes_used": len(scenes),
    "mean_area_a_c": avg("area_a_c"),
    "mean_area_r_c": avg("area_r_c"),
    "mean_difference_c": avg("difference_c"),
    "ndvi_area_a": mean_ndvi(area_a),
    "ndvi_area_r": mean_ndvi(area_r),
    "area_a_bounds": AREA_A,
    "area_r_bounds": AREA_R,
    "scenes": scenes,
}

os.makedirs("data", exist_ok=True)
with open("data/latest.json", "w") as f:
    json.dump(out, f, indent=2)
print("Wrote data/latest.json with", len(scenes), "scenes")
