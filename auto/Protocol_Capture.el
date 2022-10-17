(TeX-add-style-hook
 "Protocol_Capture"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "10pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("geometry" "margin=1in" "includefoot" "nohead") ("xcolor" "svgnames")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "fontenc"
    "amsmath"
    "amssymb"
    "amsfonts"
    "epstopdf"
    "sectsty"
    "fancyhdr"
    "graphicx"
    "booktabs"
    "makecell"
    "fancyref"
    "hyperref"
    "siunitx"
    "physics"
    "listings"
    "romannum"
    "multirow"
    "longtable"
    "geometry"
    "xcolor"
    "microtype"
    "fontspec"
    "unicode-math")
   (LaTeX-add-labels
    "sec:msa"
    "tab:access_dates"
    "tab:designed-positions"
    "sec:relax"
    "tab:epitopes")
   (LaTeX-add-xcolor-definecolors
    "Darkgray"))
 :latex)

