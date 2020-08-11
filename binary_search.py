#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 14:54:19 2020

@author: kun-linho
"""

def binarySearch (arr, left, right, search_target): 
  
    # Check base case 
    if right >= left: 
  
        mid = int(left + (right - left)/2)
  
        # If element is present at the middle itself 
        if arr[mid] == search_target: 
            return mid 
          
        # If element is smaller than mid, then it can only 
        # be present in left subarray 
        elif arr[mid] > search_target: 
            return binarySearch(arr, left, mid-1, search_target) 
  
        # Else the element can only be present in right subarray 
        else: 
            return binarySearch(arr, mid+1, right, search_target) 
  
    else: 
        # Element is not present in the array 
        return -1
    
    
arr = [ 2, 3, 4, 10, 40 ] 
x = 10
  
# Function call 
result = binarySearch(arr, 0, len(arr)-1, x) 

if result != -1: 
    print ("Element is present at index %d" % result )
else: 
    print ("Element is not present in array")