from math import pi, sqrt, sin, cos, acos
import copy
import pandas as pd

def convert(theta, R, origin_point):
    x = R*cos(theta)
    y = R*sin(theta)
    x_related = x + origin_point[0]
    y_related = y + origin_point[1]
    return x_related, y_related

def compute_convex(p1, p2):
    x1 = p1[0]
    y1 = p1[1]
    x2 = p2[0]
    y2 = p2[1]
    if x1 == x2:
        x2 -= 0.1
    k = (y2-y1)/(x2-x1)
    l = sqrt((x2-x1)**2 + (y2-y1)**2)
    d_l = l * sin(15/180*pi)
    if k==0:
        d_x = sqrt(d_l**2)
        d_y = 0
    else:
        d_x = sqrt(d_l**2/(1+1/k**2))
        d_y = -(1/k)*d_x
    
    direct = int(x1+x2)%2

    x3 = (x1+x2)/2 + (-1)**direct*d_x
    y3 = (y1+y2)/2 + (-1)**direct*d_y
    return x3, y3

def unit_angle_perimeter(R, rnode):
    min_theta = pi-2*(acos(2*rnode/R))
    k = int(2*pi/min_theta)
    theta = 2*pi/k
    return theta, k

def nlayer(max_width, rnode):
    if max_width < rnode:
        print('must have at least one node')
        return 0
    n = int((max_width/2 - rnode)/(4*rnode)) + 1
    return n

def estimate_nlayer(nnode):
    rnode = 1
    if nnode == 1:
        return 1
    nlayer = 1
    n = 1
    R = 1
    while n<nnode:
        R += 4*rnode
        theta, layer_nnode = unit_angle_perimeter(R, rnode)
        n += layer_nnode
        nlayer += 1
    return  nlayer

def min_width(rnode, nlayer):
    return ((nlayer-1)*4*rnode + rnode)*2

# return width, polar df
def assign_pos(nnode, rnode, margin=0):
    polar_df = pd.DataFrame(columns=['theta', 'r'])
    nlayer = estimate_nlayer(nnode)
    width = min_width(rnode, nlayer) + margin*2
    origin = (width/2, width/2)
    unit_angle, layer_nnode = unit_angle_perimeter((nlayer-1)*4*rnode+rnode, rnode)
    rotate_angle = unit_angle/2
    polar_df.loc[0, ] = [0, 0]
    n = 1
    for i in range(1, nlayer):
        start_angle = i * rotate_angle
        R = 4*rnode*i
        unit_angle, layer_nnode = unit_angle_perimeter(R, rnode)
        for j in range(layer_nnode):
            if n >= nnode:
                break
            theta = start_angle + unit_angle*j
            polar_df.loc[n, ] = [theta, R]
            n += 1
    return width, origin, polar_df

def polar_squar_distance(theta1,r1, theta2, r2):
    d = r1**2 + r2**2 - 2*r1*r2*cos(theta1-theta2)
    return d

def sort_by_d(source_theta, source_r, polar_df):
    polar_sorted = copy.deepcopy(polar_df)
    polar_sorted['d'] = None
    for i in polar_df.index:
        theta = polar_df.loc[i, 'theta']
        r = polar_df.loc[i, 'r']
        polar_d = polar_squar_distance(source_theta, source_r, theta, r)
        polar_df.loc[i, 'd'] = polar_d
    sorted_df = polar_df.sort_values(by=['d', 'theta', 'r'])
    return sorted_df

def ring_layout(color_dict, polar_df):
    polar_sorted = polar_df.sort_values(by=['r', 'theta'], ascending=True)
    polar_sorted['color'] = None
    order_color = sorted(color_dict, key=color_dict.get, reverse=False)
    start = 0
    for color in order_color:
        num = color_dict[color]
        polar_sorted.iloc[start: start+num, 2] = color
        start = start + num
    return polar_sorted

def sector_layout(color_dict, polar_df):
    polar_sorted = polar_df.sort_values(by=['theta', 'r'], ascending=True)
    polar_sorted['color'] = None
    order_color = sorted(color_dict, key=color_dict.get, reverse=True)
    start = 0
    for color in order_color:
        num = color_dict[color]
        polar_sorted.iloc[start: start+num, 2] = color
        start = start + num
    return polar_sorted

def scale_layout(color_dict, polar_df, width):
    r = width/4
    seeds = []
    seed = (0, r)
    order_color = sorted(color_dict, key=color_dict.get, reverse=True)
    polar_df['color'] = None
    remaining = copy.deepcopy(polar_df)
    for color in order_color:
        seeds.append(seed)
        num = color_dict[color]
        if num == 0:
            continue
        polar_sorted = sort_by_d(seed[0], seed[1], remaining) 
        selected = polar_sorted.iloc[:num, ]
        for idx in selected.index:
            polar_df.loc[idx, 'color'] = color
        remaining = polar_sorted.iloc[num:, ]
        max_theta = max(selected['theta'])
        seed = (max_theta, r)
    return polar_df, seeds
