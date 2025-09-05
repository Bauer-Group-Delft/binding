# -*- coding: utf-8 -*-
import re

def read_data(filename, runs):
    f = open(f"{filename}_data.txt", "r")
    data = []
    
    while True:
        line = f.readline()
        if not line:
            # End of file
            break
        
        # Obtain keys
        keys = re.split(', *', line.replace('\n', ''))
        
        # Begin next dataset
        data.append({})
        for key in keys:
            data[-1][key] = []
            
        # Store data
        for i in range(runs):
            line = f.readline()
            content = re.split(', *', line.replace('\n', ''))
            for k,key in enumerate(keys):
                data[-1][key].append(content[k])
 
    f.close()
    return data
