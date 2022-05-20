(TeX-add-style-hook
 "misclassifications"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("babel" "english")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "babel"
    "booktabs"
    "caption"
    "makecell"
    "multirow"
    "amsmath"))
 :latex)

