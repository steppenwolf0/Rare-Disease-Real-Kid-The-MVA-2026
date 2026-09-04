#!/usr/bin/env python3

import requests
import pandas as pd

# Queries to run
QUERIES = [
    "readthrough therapy",
    "premature stop codon",
    "nonsense mutation",
    "Mosaic Variegated Aneuploidy"
]

API_V2_BASE = "https://clinicaltrials.gov/api/v2/studies"

all_rows = []

for query in QUERIES:
    print(f"Searching: {query}")

    params = {
        "query.term": query,
        "pageSize": 1000,
        "format": "json",
    }

    r = requests.get(API_V2_BASE, params=params, timeout=30)
    r.raise_for_status()

    studies = r.json().get("studies", [])

    for st in studies:
        ident = st.get("protocolSection", {}).get("identificationModule", {})
        status = st.get("protocolSection", {}).get("statusModule", {})
        design = st.get("protocolSection", {}).get("designModule", {})

        nct_id = ident.get("nctId", "")
        title = ident.get("briefTitle", "")

        all_rows.append({
            "query": query,
            "nct_id": nct_id,
            "title": title,
            "status": status.get("overallStatus", ""),
            "phase": "; ".join(design.get("phases", [])),
            "url": f"https://clinicaltrials.gov/study/{nct_id}" if nct_id else "",
        })

df = pd.DataFrame(all_rows)
df.to_csv("clinicaltrials_results.csv", index=False)

print(f"Found {len(df)} studies")