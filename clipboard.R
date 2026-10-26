clipboard <- htmltools::tagList(
  xaringanExtra::use_clipboard(
    button_text = "<i class=\"fa-solid fa-clipboard\" style=\"color: #00008B\"></i>",
    success_text = "<i class=\"fa fa-check\" style=\"color: #90BE6D\"></i>",
    error_text = "<i class=\"fa fa-times-circle\" style=\"color: #F94144\"></i>"
  ),
  rmarkdown::html_dependency_font_awesome()
)