import pandas as pd
import numpy as np
import numexpr
def string_search(data, column_name,item, reverse = False):
    if reverse == False:
        return data[data[column_name].to_numpy() == item]
    else:
        return data[data[column_name].to_numpy() != item]
def num_search(data, column_name,number, direction, step = None,inclusion = False):
    x = data[column_name].values
    if direction == ">":
        if inclusion == False:
            return(data[numexpr.evaluate('(x > number)')])
        else:
            return(data[numexpr.evaluate('(x >= number)')])
    elif direction == '<':
        if inclusion == False:
            return(data[numexpr.evaluate('(x < number)')])
        else:
            return(data[numexpr.evaluate('(x <= number)')])

    elif direction == '==':
        return(data[numexpr.evaluate('(x == number)')])
    elif direction =='between' and step != None:
        if inclusion == False:
            temp = data[numexpr.evaluate('(x > number-step)')]
            x = temp[column_name].values
            return (temp[numexpr.evaluate('(x < number+step)')])
        else:
            temp = data[numexpr.evaluate('(x >= number-step)')]
            x = temp[column_name].values
            return (temp[numexpr.evaluate('(x <= number+step)')])
    else:
        print('the wrong method is passed')



