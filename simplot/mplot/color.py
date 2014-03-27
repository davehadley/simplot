import colorsys

def brighten(rgb, mag=None):
    return _change_brightness(rgb, +1, mag)

def darken(rgb, mag=None):
    return _change_brightness(rgb, -1, mag)

def _change_brightness(rgb, sign, mag=None):
    hsv = colorsys.rgb_to_hsv(*rgb)
    v = hsv[2]
    if mag is None:
        mag = 0.5
    #constrain magnitude to 0-1
    mag = max(0.0, mag)
    mag = min(1.0, mag)
    #calcualte increment
    incr = 0.0
    if sign > 0:
        incr = (1.0 - v) * mag
    elif sign < 0:
        incr = (v - 0.0) * -mag
    newv = v + incr
    #ensure result is on the line 0-1
    newv = min(1.0, newv)
    newv = max(0.0, newv)
    rgb = colorsys.hsv_to_rgb(hsv[0], hsv[1], newv) 
    return rgb