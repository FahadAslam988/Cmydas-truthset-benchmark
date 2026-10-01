"""Shared style for every Paper 1 figure (one palette, one font, one set of sizes)."""
import matplotlib
from matplotlib import font_manager
# Journal-standard scheme (Nature figure guide -> Wong 2011, Nature Methods 8:441; Okabe-Ito palette):
# categorical = Okabe-Ito exactly; ordered values = one-hue sequential (ColorBrewer Blues); diagrams = greys + one accent.
NAVY, TEAL, MINT, PALE, GREY_TXT = '#08519C', '#3182BD', '#6BAED6', '#F2F4F7', '#3A3A3A'   # Blues family
CALLERS = ['BCF', 'Deepvariant', 'Freebayes', 'Gatk', 'Varscan']
LABEL = {'BCF': 'BCFtools', 'Deepvariant': 'DeepVariant', 'Freebayes': 'FreeBayes', 'Gatk': 'GATK', 'Varscan': 'VarScan'}
COL = {'BCF': '#E69F00', 'Deepvariant': '#0072B2', 'Freebayes': '#009E73', 'Gatk': '#CC79A7', 'Varscan': '#D55E00'}  # Okabe-Ito
MK = {'BCF': 'o', 'Deepvariant': 's', 'Freebayes': '^', 'Gatk': 'D', 'Varscan': 'v'}
RANK_SHADE = {1: '#08519C', 2: '#3182BD', 3: '#6BAED6', 4: '#BDD7E7', 5: '#EFF3FF'}  # ColorBrewer Blues
FS, FS_SMALL, FS_HEAD, FS_LETTER = 7, 6, 7, 9
def setup():
    fam = 'Arial' if any('Arial' in f.name for f in font_manager.fontManager.ttflist) else 'Liberation Sans'
    matplotlib.rcParams.update({'font.family': fam, 'font.size': FS, 'axes.linewidth': 0.8, 'pdf.fonttype': 42,
        'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'axes.labelcolor': '#1F1F1F', 'text.color': '#1F1F1F'})
def clean(ax, grid=True):
    ax.spines[['top', 'right']].set_visible(False)
    if grid: ax.grid(axis='y', color='0.92', lw=0.6); ax.set_axisbelow(True)
def letter(ax, x, y, s):
    ax.text(x, y, s, transform=ax.transAxes, fontweight='bold', fontsize=FS_LETTER, va='bottom')
def head(ax, x, y, s, **kw):
    ax.text(x, y, s, ha='center', va='bottom', fontsize=FS_HEAD, fontweight='bold', color=NAVY, **kw)
def save(fig, out, name):
    for ext, kw in (('pdf', {}), ('eps', {}), ('tiff', {'dpi': 300, 'pil_kwargs': {'compression': 'tiff_lzw'}}), ('png', {'dpi': 300})):
        fig.savefig(f'{out}/{name}.{ext}', bbox_inches='tight', facecolor='white', **kw)
