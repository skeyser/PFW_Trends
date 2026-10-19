#!/bin/bash
#Daymet Constants
years="1980,1981,1982,1983,1984,1985,1986,1987,1988,1989"
years+=",1990,1991,1992,1993,1994,1995,1996,1997,1998,1999"
years+=",2000,2001,2002,2003,2004,2005,2006,2007,2008,2009"
years+=',2010,2011,2012,2013,2014,2015,2016,2017,2018,2019'
years+=',2020,2021,2022,2023,2024,2025'
vars="dayl,prcp,srad,swe,tmax,tmin,vp"
# read in first command line argument (should be latlon.txt)
inF=$1

# check if first command line argument exists
if [ -z $inF ]
then
  echo Error: No input file provided
  echo Usage: ./run.sh inF
  exit 1
fi

inFs=()
lats=()
lons=()
count=0

while IFS= read -r p; do
  count=$((count + 1))

  if (( count % 1000 == 0 )); then
    echo "Read $count input lines..." >&2
  fi

  # Remove whitespace
  p=${p//[[:space:]]/}

  # Ignore blank lines
  [[ -z "$p" ]] && continue

  # Optional Years:... line
  if [[ "$p" =~ ^[Yy]ears: ]]; then
    years=${p#*:}
    continue
  fi

  # Optional Variables:... line
  if [[ "$p" =~ ^[Vv]ariables: ]]; then
    vars=${p#*:}
    continue
  fi

  # Read filename, latitude, longitude
  IFS=',' read -r outfile lat lon <<< "$p"

  if [[ -n "$lon" ]]; then
    inFs+=("$outfile")
    lats+=("$lat")
    lons+=("$lon")
  else
    # Preserve original two-column behavior: latitude,longitude
    inFs+=("NULL")
    lats+=("$outfile")
    lons+=("$lat")
  fi
done < "$inF"

echo "Finished reading $count input lines. Starting downloads..." >&2

if [ -x "$(command -v curl)" ]; then
  downloadwName() { curl -J -O --silent $1; }
  download() { curl --silent $1 > $2; }
else
  if [ -x "$(command -v wget)" ]; then
    downloadwName() { wget --content-disposition -q $1; }
    download() { wget -q $1 -o $2; }
  else
    echo "Error: You must have either curl or wget installed to run this program."
    exit 1
  fi
fi

max_jobs=${MAX_JOBS:-4}

export years vars
export -f download downloadwName

for key in "${!lats[@]}"; do
  printf '%s\0%s\0%s\0' \
    "${inFs[$key]}" \
    "${lats[$key]}" \
    "${lons[$key]}"
done |
xargs -0 -n 3 -P "$max_jobs" bash -c '
  outfile=$1
  lat=$2
  lon=$3

  url="https://daymet.ornl.gov/single-pixel/api/data?lat=${lat}&lon=${lon}&measuredParams=${vars}&year=${years}"

  echo "Processing $url"

  if [[ "$outfile" != "NULL" ]]; then
    download "$url" "$outfile"
  else
    downloadwName "$url"
  fi
' _

echo ""
