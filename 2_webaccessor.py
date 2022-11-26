# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 18:21:51 2022

@author: danie
"""
import time
from os.path import exists
from pathlib import Path
from sys import exit
from selenium import webdriver
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import copy
import os
import shutil


SCOPnames = {}
SCOPnames1 = []
file1 = open('Human-Readable SCOP_entry.txt', 'r')
Lines = file1.readlines()

for line in Lines:
    if line[1:line.find(' ')].isnumeric():
        SCOPnames1.append(line[1:line.find(' ')])
    elif line[line.rfind(' '):-1].isnumeric():
        SCOPnames1.append(line[line.rfind(' '):-1])
SCOPnames = set(SCOPnames1)
# print(SCOPnames)
PDBnames = []
try:
    infile = open('PDBnames.txt','r')
    temp = infile.read().splitlines()
    infile.close()
    PDBnames = temp
except:
    
    PDBnames = []
    for name in SCOPnames:
        
        driver = webdriver.Firefox(executable_path=r'C:\Users\danie\Downloads\geckodriver-v0.30.0-win64\geckodriver.exe')
        driver.get("https://scop.mrc-lmb.cam.ac.uk/")
        
        
        inputElement = driver.find_element_by_xpath("//*[@id=\"searchbox\"]")
        inputElement.send_keys('%s'%name)
        inputElement.send_keys(Keys.ENTER)
        
        
        timeout = 15
        try:
            element_present = EC.presence_of_element_located((By.XPATH, "/html/body/div/div/div[2]/div[2]/div[2]/div/div/div/a"))
            WebDriverWait(driver, timeout).until(element_present)
        except TimeoutException:
            print("Timed out waiting for page to load")
    
        PDBname = driver.find_element_by_xpath("/html/body/div/div/div[2]/div[2]/div[2]/div/div/div/a").text[0:4]
        PDBnames.append(PDBname)
        driver.close()
    outfile = open('PDBnames.txt','w')
    for name in PDBnames:
        outfile.write(name)
        outfile.write('\n')
    outfile.close()

print(PDBnames)



# try:
#     infile = open('PDBnames.txt','r')
#     temp = infile.read().splitlines()
#     infile.close()
#     PDBnames = temp
# except:

    

currentDir = os.getcwd()
os.chdir('C:/Users/danie/Downloads')

for name in PDBnames:
    try:
        exists('C:/Users/danie/Downloads/%s'%(name.lower()+'.pdb'))
    except:
        print('The file %s did not exist in downloads folder'%(name.lower()+'.pdb'))

        for name in PDBnames:
            driver = webdriver.Firefox(executable_path=r'C:\Users\danie\Downloads\geckodriver-v0.30.0-win64\geckodriver.exe')
            driver.get("https://www.rcsb.org/")
            
            
            inputElement = driver.find_element_by_xpath("//*[@id=\"search-bar-input-text\"]")
            inputElement.send_keys(name)
            time.sleep(0.5)
            clickies = driver.find_element_by_xpath('/html/body/div/div[2]/div/div/div[1]/div[2]/div/div[2]/table/tbody/tr[1]/td[3]/div/span').click()
            time.sleep(0.5)
            timeout = 15
            try:
                element_present = EC.presence_of_element_located((By.XPATH, "//*[@id=\"dropdownMenuDownloadFiles\"]"))
                WebDriverWait(driver, timeout).until(element_present)
            except TimeoutException:
                print("Timed out waiting for page to load")
            
            clickies = driver.find_element_by_xpath("//*[@id=\"dropdownMenuDownloadFiles\"]").click()
            
            time.sleep(0.1)
            clickies = driver.find_element_by_xpath("/html/body/div/div[3]/div[2]/div[1]/div/div[2]/div[1]/div[2]/ul/li[3]/a").click()
            time.sleep(0.1)

for name in PDBnames:
    driver = webdriver.Firefox(executable_path=r'C:\Users\danie\Downloads\geckodriver-v0.30.0-win64\geckodriver.exe')
    driver.get("http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py")
    
    
    inputElement = driver.find_element("xpath", "/html/body/div/div/form/table[1]/tbody/tr[2]/td[2]/input")
    
    path_to_download_folder = str(os.path.join(Path.home(), "Downloads"))
    
    
    inputElement.send_keys(path_to_download_folder + "\%s.pdb"%name)
    inputElement = driver.find_element("xpath",'/html/body/div/div/form/table[1]/tbody/tr[3]/td[2]/input')
    time.sleep(0.1)
    inputElement.send_keys(Keys.ENTER)
    time.sleep(5)
    
    
    file = driver.find_element("xpath",'/html/body').text
    
    with open('%s.txt'%name, 'w', encoding='utf-8') as outfile:
        print(file, file=outfile)
    time.sleep(0.1)
    os.replace("C:/Users/danie/Downloads/%s.txt"%name, "C:/Users/danie/Desktop/capstone/v3/%s.txt"%name)
    driver.close()

exit()



# /html/body/div/div/form/table[1]/tbody/tr[3]/td[2]/input