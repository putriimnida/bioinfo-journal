---
title: "bulkRNA-seq"
author: "Putri Ramadani"
date: "2026-11-23"
output: html_document
---


# Verifying FASTQ files against their MD5 checksums
An MD5 checksum is a 32 character hexadecimal number that is computed on a file. 
If two files have the same MD5 checksum value, it is highly probable that the two files are the same.
It is generally used to check data integrity.
For example: d41d8cd98f00b204e9800998ecf1000e.
Change even a single byte of the file and the checksum comes out completely different. 
That makes it useful for checking that a file arrived intact.

## for Mac users, use md5sum filename
```bash
# go to the data_release directory where fastq is at

cat MD5.txt        # prints the whole file
less MD5.txt       # scroll through it (best for big files)
md5sum MD5.txt     # calculate the checksum
md5sum -c MD5.txt   
```
With -c, each file is reported as OK or FAILED.
One caveat: MD5 is fine for catching accidental corruption. 
It's no longer considered secure against deliberate tampering, because people can create two different files with the same MD5. 
For security-sensitive checks, use SHA-256 (sha256sum).

## for Windows users
- Double click "MD5-win.exe", click OK on the popped window
- Click "Verify sums" and choose "MD5.txt" in data_release folder, then click OK
- Wait until the task is done
Red indicates error
Yellow indicates processing
Green indicates checking pass

Source: Novogene MD5 manual




# FastQC on the raw reads
get a .html report for reads quality, possible contamination etc
```bash
module load fastqc
mkdir -p fastqc_out
fastqc -o fastqc_out \
  /cluster/projects/xxx/01.RawData/CAF2_10/*.fq.gz \ #lowreads
  /cluster/projects/xxx/01.RawData/HTSL_3/*.fq.gz #normal sample with high reads
```




