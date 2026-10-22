import sys, requests, json, time

vcf_in = sys.argv[1]
csv_out = sys.argv[2]

queries = []
with open(vcf_in, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        cols = line.strip().split('\t')
        
        # FIX: Strip 'chr' prefix because Ensembl API expects '6' instead of 'chr6'
        chrom = cols[0].replace('chr', '')
        start = int(cols[1])
        ref = cols[3]
        alt = cols[4].split(',')[0]
        end = start + len(ref) - 1
        
        queries.append(f"{chrom} {start} {end} {alt} 1")

if not queries:
    print("No PASS variants found. Exiting.")
    sys.exit(0)

print(f"[ANNOTATE] Sending {len(queries)} variants to Ensembl VEP REST API...")
server = "https://rest.ensembl.org"
ext = "/vep/human/region"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

results = []
for i in range(0, len(queries), 200):
    payload = {"variants": queries[i:i+200]}
    success = False
    
    for attempt in range(5):
        r = requests.post(server+ext, headers=headers, data=json.dumps(payload))
        if r.ok:
            results.extend(r.json())
            success = True
            break
        elif r.status_code >= 500 or r.status_code == 429:
            print(f"  [WARN] API Overloaded (HTTP {r.status_code}). Retrying in 10s... (Attempt {attempt+1}/5)")
            time.sleep(10)
        else:
            print(f"  [ERROR] API Error {r.status_code}: {r.text}")
            break
            
    if not success:
        print("  [FATAL] Ensembl API is completely down.")
        sys.exit(1)
        
    time.sleep(0.5)

with open(csv_out, 'w') as out:
    out.write("Location,Allele,Gene,Consequence,Impact\n")
    for res in results:
        loc = f"chr{res.get('seq_region_name', '')}:{res.get('start', '')}"
        allele = res.get('allele_string', '')
        
        tv_list = res.get('transcript_consequences', [{}])
        tv = tv_list[0] if tv_list else {}
        
        gene = tv.get('gene_symbol', tv.get('gene_id', 'Unknown'))
        conseq = tv.get('consequence_terms', ['Unknown'])[0]
        impact = tv.get('impact', 'Unknown')
        
        out.write(f"{loc},{allele},{gene},{conseq},{impact}\n")
        
print(f"  -> Successfully saved functional annotations to {csv_out}")
