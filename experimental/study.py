import requests
from bs4 import BeautifulSoup
from pprint import pprint
import sys
import csv

def __main__():
    infile = "TSX_cleaned.txt"
    symbols = get_symbols(infile)
    #symbols = ["ILMN", "JD"]
    epss = dict()
    for symbol in symbols:
        eps = get_finance(symbol)
        epss[symbol] = eps
    pprint(epss)
    for symbol in epss:
        items = [i for i in epss[symbol]]
        item_list = items[0] + items[1]
        print symbol, "\t".join(item_list)

def get_symbols(infile):
    """ Dictionary holds all files: patient -> status -> file_identifier -> file_path  """
    symbols = []
    with open(infile, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            symbols.append(line['Symbol'])
    print symbols
    return symbols


def get_finance(symbol):
    url = "http://www.google.com/finance?q=%s&fstype=ii" % symbol
    #r = requests.get("http://www.google.com/finance?q=%s&fstype=ii") % symbol
    r = requests.get(url)
    soup = BeautifulSoup(r.content)
    tables = soup.find_all('table')
    eps = dict()
    eps_vals = []
    for table in tables:
        results = {}
        th = table.find_all('tr')#,text=['menu1','menu2'])
        for x in th:
            #print x
            results_li = []
            li = x.find_all('td')[1:]
            for y in li:
                #print y.next
                results_li.append(y.text)
            results[x.text.split('\n')[1]] = results_li
        #pprint(results)
        for key in results:
            if ('Diluted Normalized EPS' in key):
                eps_vals.append(results[key])

    return eps_vals
    
if __name__ == '__main__':
    __main__()

sys.exit()

# below code for debug
#i = 1
#for table in tables:
#    i=i+1
#    print "-------------table %s---------------------" % i
#    #print table
#    trs = table.find_all('tr')
#    for tr in trs:
#       tds = tr.find_all('td')
#       for td in tds:
#           print td.text            
#           #print "--------"
#
#
#tables = soup.find_all(lambda tag: tag.name=='table' and tag.has_attr('id') and tag['id']=="fs-table")

#for table in tables:
#    rows = table.find_all(lambda tag: tag.name=='tr')
#    print rows[-1]
