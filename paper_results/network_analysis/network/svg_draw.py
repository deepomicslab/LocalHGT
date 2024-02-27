# use python code to draw svg
def canvas(width, height, elements, style=""):
    svg_text = '\
<svg version="1.1"\n \
        baseProfile="full"\n \
        width="{}" height="{}"\n \
        xmlns="http://www.w3.org/2000/svg">\n \
    <style>\n \
{} \n\
    </style>\n \
    {}\n</svg>'.format(width, height, style, elements)
    return svg_text

def draw_circle(x, y, r="10", stroke="black", fill="none", swidth="1", cname="circle", id="", style=False):
    if not style:
        if id == "":
            circle_tmp = '<circle class="{}" cx="{}" cy="{}"/>\n'
            circle = circle_tmp.format(cname, x, y)
        else:
            circle_tmp = '<circle id="{}" class="{}" cx="{}" cy="{}"/>\n'
            circle = circle_tmp.format(id, cname, x, y)
        return circle
    if id == "":
        circle_tmp = '<circle class="{}" cx="{}" cy="{}" r="{}" stroke="{}" stroke-width="{}"  fill="{}" />\n'
        circle = circle_tmp.format(cname, x, y, r, stroke, swidth, fill)
    else:
        circle_tmp = '<circle id="{}" class="{}" cx="{}" cy="{}" r="{}" stroke="{}" stroke-width="{}"  fill="{}" />\n'
        circle = circle_tmp.format(id, cname, x, y, r, stroke, swidth, fill)
    return circle

def draw_curve(x1, y1, x2, y2, x3, y3, stroke="black", fill="none", swidth="1", cname="curve", id="", style=False):
    if not style:
        if id == "":
            path_tmp = '<path class="{}" d="M {} {} Q {} {} {} {}" />\n'
            curve = path_tmp.format(cname, x1, y1, x2, y2, x3, y3)
        else:
            path_tmp = '<path id="{}" class="{}" d="M {} {} Q {} {} {} {}" />\n'
            curve = path_tmp.format(id, cname, x1, y1, x2, y2, x3, y3)

        return curve
    if id == "":
        path_tmp = '<path class="{}" d="M {} {} Q {} {} {} {}" stroke="{}" stroke-width="{}" fill="{}" />\n'
        curve = path_tmp.format(cname, x1, y1, x2, y2, x3, y3, stroke, swidth, fill)
    else:
        path_tmp = '<path id="{}" class="{}" d="M {} {} Q {} {} {} {}" stroke="{}" stroke-width="{}" fill="{}" />\n'
        curve = path_tmp.format(id, cname, x1, y1, x2, y2, x3, y3, stroke, swidth, fill)
    return curve

def draw_text(content, x, y, anchor="start", cname="text", id=""):
    if id == "":
        text_tmp = '<text text-anchor="{}" x="{}" y="{}" class="{}">{}</text>\n'
        text = text_tmp.format(anchor, x, y, cname, content)
    else:
        text_tmp = '<text id="{}" text-anchor="{}" x="{}" y="{}" class="{}">{}</text>\n'
        text = text_tmp.format(id, anchor, x, y, cname, content)
    return text
    