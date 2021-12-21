from matplotlib.text import TextPath
from matplotlib.font_manager import FontProperties

fp = FontProperties(family="Arial", weight="bold")
globscale = 1.35

color_schema = {
    'A':'k','C':'g','D':'r','E':'r','F':'k','G':'g','H':'b','I':'k','K':'b','L':'k',
    'M':'k','N':'pink','P':'k','Q':'pink','R':'b','S':'b','T':'b','V':'k','W':'k','Y':'g'}

letter_schema = {
    'A':TextPath((-0.35,0),'A',size=1,prop=fp),
    'C':TextPath((-0.35,0),'C',size=1,prop=fp),
    'D':TextPath((-0.35,0),'D',size=1,prop=fp),
    'E':TextPath((-0.35,0),'E',size=1,prop=fp),
    'F':TextPath((-0.30,0),'F',size=1,prop=fp),
    'G':TextPath((-0.38,0),'G',size=1,prop=fp),
    'H':TextPath((-0.35,0),'H',size=1,prop=fp),
    'I':TextPath((-0.14,0),'I',size=1,prop=fp),
    'K':TextPath((-0.38,0),'K',size=1,prop=fp),
    'L':TextPath((-0.33,0),'L',size=1,prop=fp),
    'M':TextPath((-0.42,0),'M',size=1,prop=fp),
    'N':TextPath((-0.36,0),'N',size=1,prop=fp),
    'P':TextPath((-0.34,0),'P',size=1,prop=fp),
    'Q':TextPath((-0.39,0),'Q',size=1,prop=fp),
    'R':TextPath((-0.37,0),'R',size=1,prop=fp),
    'S':TextPath((-0.33,0),'S',size=1,prop=fp),
    'T':TextPath((-0.31,0),'T',size=1,prop=fp),
    'V':TextPath((-0.33,0),'V',size=1,prop=fp),
    'W':TextPath((-0.47,0),'W',size=1,prop=fp),
    'Y':TextPath((-0.33,0),'Y',size=1,prop=fp),
}