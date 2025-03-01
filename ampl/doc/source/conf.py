# -*- coding: utf-8 -*-
#
import sys
import os
import re
import datetime

# -- General configuration -----------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "3.2.0"

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.mathjax",
    "sphinxcontrib.googleanalytics",
    "sphinx_sitemap",
]

# Configure Breathe.
# When building with CMake, the path to doxyxml is passed via the command line.
breathe_projects = {"mp": "doxyxml"}
breathe_default_project = "mp"
breathe_domain_by_extension = {"h": "cpp"}

highlight_language = "c++"
primary_domain = "cpp"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix of source filenames.
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "amplgsl"
copyright = "2016-{}, AMPL Optimization Inc".format(datetime.date.today().year)

# The short X.Y version.
version = "1.0"

# The full version, including alpha/beta/rc tags.
release = "20211111"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "virtualenv"]

# The reST default role (used for this markup: `text`) to use for all documents.
default_role = "cpp:any"

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "ampl_sphinx_theme"
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/ampl/gsl",
            "icon": "fab fa-github",
        },
    ],
    "logo_text": "GSL",
}
html_context = {"default_mode": "light"}
googleanalytics_id = "G-0D29096DY8"


html_baseurl = "https://gsl.ampl.com"

html_extra_path = ["_html"]

sitemap_filename = "sphinx-sitemap.xml"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
# html_css_files = [
#     'css/custom.css',
# ]
# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = "amplgsl {} documentation".format(release)
html_title = "amplgsl documentation"

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "amplgsl"

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = False

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "https://raw.githubusercontent.com/ampl/ampl.github.io/master/themes/static/ampl-navbar-logo.png"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = None
html_favicon = "https://raw.githubusercontent.com/ampl/ampl.github.io/master/themes/static/ampl-favicon.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = "AMPLdoc"


# -- Options for LaTeX output --------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
    (
        "index",
        "amplgsl.tex",
        "AMPLGSL Documentation",
        "AMPL Optimization, Inc.",
        "manual",
    ),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
# latex_use_parts = False

# If true, show page references after internal links.
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
# latex_domain_indices = True


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [("index", "ampl", "AMPLGSL Documentation", ["AMPL Optimization, Inc."], 1)]

# If true, show URL addresses after external links.
# man_show_urls = False


# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        "index",
        "AMPL",
        "AMPLGSL Documentation",
        "AMPL Optimization, Inc.",
        "AMPL",
        "Documentation for AMPL GSL bindings.",
        "Miscellaneous",
    ),
]

# Documents to append as an appendix to all manuals.
# texinfo_appendices = []

# If false, no module index is generated.
# texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
# texinfo_show_urls = 'footnote'


def builder_inited_handler(app):
    import os
    import sys

    sys.path.append(os.path.dirname(os.path.realpath(__file__)))
    from support import extractdocs

    me = os.path.dirname(os.path.realpath(__file__))
    src = os.path.join(me, "../../src/amplgsl.cc")
    dest = os.path.join(me, "ref")
    extractdocs.extract_docs(src, dest)


def setup(app):
    app.connect("builder-inited", builder_inited_handler)
