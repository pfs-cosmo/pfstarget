#!/usr/bin/env python
# Based on
# https://hsc-gitlab.mtk.nao.ac.jp/ssp-software/data-access-tools/-/blob/master/dr3/catalogQuery/hscSspQuery3.py
import os
import sys
import csv
import json
import time
import astropy.io.fits as pyfits
import getpass
import argparse
import urllib.request, urllib.error, urllib.parse
import astropy.io.ascii as ascii

version =   20190924.1
args    =   None
doDownload  =   True
doUnzip =   True
diffver =   '-colorterm'

def chunkNList(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--user', '-u', required=True,
                        help='specify your STARS account')
    parser.add_argument('--delete-job', '-D', action='store_true',
                        help='delete job after your downloading')
    parser.add_argument('--format', '-f', dest='out_format', default='fits',
                        choices=['csv', 'csv.gz', 'sqlite3', 'fits'],
                        help='specify output format')
    parser.add_argument('--nomail', '-M', action='store_true',
                        help='suppress email notice')
    parser.add_argument('--release_year', '-r',
                        choices='s15 s16 s17 s18 s19 s20 s21 s23'.split(),
                        required=True,
                        help='the release year')
    parser.add_argument('--password-env', default='HSC_SSP_CAS_PASSWORD',
                        help='environment variable for STARS')
    parser.add_argument('--preview', '-p', action='store_true',
                        help='quick mode (short timeout)')
    parser.add_argument('--skip-syntax-check', '-S', action='store_true',
                        help='skip syntax check')
    parser.add_argument('--api-url',
                        default='https://hscdata.mtk.nao.ac.jp/datasearch/api/catalog_jobs/',
                        help='for developers')
    parser.add_argument('sql-file', type=argparse.FileType('r'),
                        help='SQL file')

    global args,release_version,prefix,prefix2,ngroups
    args = parser.parse_args()
    release_year   =   args.release_year
    prefix  =   'database/%s%s/sql' %(release_year,diffver)
    prefix2 =   'database/%s%s/tracts' %(release_year,diffver)
    if not os.path.exists(prefix):
        os.system('mkdir -p %s' %prefix)
    if not os.path.exists(prefix2):
        os.system('mkdir -p %s' %prefix2)

    if 's17' in release_year:
        release_version =   'dr2'
        ngroups =   5
        tractname='fieldTractInfoS16.csv'
    elif 's18' in release_year:
        release_version =   'dr2'
        ngroups =   10
        tractname='fieldTractInfoS18.csv'
    elif 's19' in release_year:
        release_version =   'dr3'
        ngroups =   20
        tractname=  'fieldTractInfoS19.csv'
    elif 's21' in release_year:
        release_version =   'dr4-citus'
        ngroups =   20
        tractname=  'fieldTractInfoS21.csv'
    elif 's23' in release_year:
        release_version =   'dr4'
        ngroups =   40
        tractname=  'TractInfoS23.csv'
    

    global sql
    sql         =   args.__dict__['sql-file'].read()
    tracts      =   ascii.read(tractname)['tract']
    tracts2     =   chunkNList(tracts,ngroups)
    if doDownload:
        global credential
        credential  =   {'account_name': args.user, 'password': getPassword()}
        for ig,tractL in enumerate(tracts2):
            #restart
            #if(ig<194): continue
            print('Group: %s' %ig)
            downloadTracts(ig,tractL)

    if doUnzip:
        for ig,tractL in enumerate(tracts2):
            print('unzipping group: %s' %ig)
            separateTracts(ig,tractL)
    
    return

def separateTracts(ig,tractL):
    infname     =   '%s.%s'%(ig,args.out_format)
    infname     =   os.path.join(prefix,infname)
    if not os.path.exists(infname):
        print('Does not have input file')
        return
    fitsAll     =   pyfits.getdata(infname)
    print('read %s galaxies' %len(fitsAll))

    for tract in tractL:
        outfname    =   os.path.join(prefix2,'%s.fits' %(tract))
        if os.path.exists(outfname):
            print('already have file for tract: %s' \
                    %tract)
            continue
        fits        =   fitsAll[fitsAll['tract']==tract]
        if len(fits)>10:
            pyfits.writeto(outfname,fits)
        del fits
    return

def downloadTracts(ig,tractL):
    tractStr    =   map(str,tractL)
    tname       =   "'{0}'".format("', '".join(tractStr))
    print(tname)
    job         =   None
    sqlU        =   sql.replace('{$tract}',tname)
    outfname    =   '%s.%s'%(ig,args.out_format)
    outfname    =   os.path.join(prefix,outfname)
    if os.path.exists(outfname):
        print('already have output')
        return
    print('querying data')
    job         =   submitJob(credential, sqlU, args.out_format)
    blockUntilJobFinishes(credential, job['id'])
    print('downloading data')
    fileOut     =   open(outfname,'w')
    fileBuffer  =   fileOut.buffer
    download(credential, job['id'], fileBuffer)
    if args.delete_job:
        deleteJob(credential, job['id'])
    print('closing output file')
    fileBuffer.close()
    fileOut.close()
    del fileBuffer
    del fileOut
    return

class QueryError(Exception):
    pass

def httpJsonPost(url, data):
    data['clientVersion'] = version
    postData = json.dumps(data)
    return httpPost(url, postData, {'Content-type': 'application/json'})

def httpPost(url, postData, headers):
    req = urllib.request.Request(url, postData.encode('utf-8'), headers)
    res = urllib.request.urlopen(req)
    return res

def submitJob(credential, sql, out_format):
    url = args.api_url + 'submit'
    catalog_job = {
        'sql'                     : sql,
        'out_format'              : out_format,
        'include_metainfo_to_body': True,
        'release_version'         : release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job, 'nomail': args.nomail, 'skip_syntax_check': args.skip_syntax_check}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job

def jobStatus(credential, job_id):
    url = args.api_url + 'status'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job

def jobCancel(credential, job_id):
    url = args.api_url + 'cancel'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)

def preview(credential, sql, out):
    url = args.api_url + 'preview'
    catalog_job = {
        'sql'             : sql,
        'release_version' : release_version,
    }
    postData=   {'credential': credential, 'catalog_job': catalog_job}
    res     =   httpJsonPost(url, postData)
    result  =   json.load(res)

    writer = csv.writer(out)
    # writer.writerow(result['result']['fields'])
    for row in result['result']['rows']:
        writer.writerow(row)

    if result['result']['count'] > len(result['result']['rows']):
        raise QueryError('only top %d records are displayed !' % len(result['result']['rows']))

def blockUntilJobFinishes(credential, job_id):
    max_interval = 0.5 * 60 # sec.
    interval = 1
    while True:
        time.sleep(interval)
        job = jobStatus(credential, job_id)
        if job['status'] == 'error':
            raise QueryError('query error: ' + job['error'])
        if job['status'] == 'done':
            break
        interval *= 2
        if interval > max_interval:
            interval = max_interval
    print('blocking over')
    return

def download(credential, job_id, out):
    url     =   args.api_url + 'download'
    postData=   {'credential': credential, 'id': job_id}
    res     =   httpJsonPost(url, postData)
    bufSize =   1024 * 1<<10 # 1024K
    bufLim  =   1 * 1<<10 # 1K
    while True:
        buf = res.read(bufSize)
        out.write(buf)
        if len(buf) < bufLim:
            break
    return

def deleteJob(credential, job_id):
    url = args.api_url + 'delete'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)
    return

def getPassword():
    password_from_envvar = os.environ.get(args.password_env, '')
    if password_from_envvar != '':
        return password_from_envvar
    else:
        return getpass.getpass('password? ')

if __name__ == '__main__':
    main()
