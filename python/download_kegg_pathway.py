#!/usr/bin/env python3
"""Download a KEGG pathway gene list via the KEGG REST API."""

import csv
import re
import sys
import urllib.request


def fetch_pathway(pathway_id: str) -> str:
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    print(f"Fetching {url} ...")
    with urllib.request.urlopen(url, timeout=30) as resp:
        return resp.read().decode("utf-8")


def parse_genes(text: str):
    genes = []
    in_gene_section = False
    for line in text.splitlines():
        if not line:
            continue
        if line[0] != " ":
            parts = line.split(None, 1)
            section = parts[0]
            if section == "GENE":
                in_gene_section = True
                content = parts[1] if len(parts) > 1 else ""
            else:
                in_gene_section = False
                continue
        elif in_gene_section:
            content = line.strip()
        else:
            continue
        m = re.match(r"^(\d+)\s+(\S+?);\s*(.+)$", content)
        if not m:
            continue
        entrez_id, symbol, description = m.groups()
        description = re.sub(r"\s*\[(?:KO|EC):[^\]]+\]", "", description).strip()
        genes.append((entrez_id, symbol, description))
    return genes


def main():
    pathway_id = sys.argv[1] if len(sys.argv) > 1 else "hsa04151"
    text = fetch_pathway(pathway_id)
    genes = parse_genes(text)
    print(f"Parsed {len(genes)} genes from {pathway_id}")

    csv_path = f"kegg_{pathway_id}_genes.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["entrez_id", "symbol", "description"])
        w.writerows(genes)
    print(f"Wrote {csv_path}")

    txt_path = f"kegg_{pathway_id}_symbols.txt"
    with open(txt_path, "w") as f:
        for _, symbol, _ in genes:
            f.write(symbol + "\n")
    print(f"Wrote {txt_path}")


if __name__ == "__main__":
    main()