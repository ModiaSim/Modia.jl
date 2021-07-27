using Documenter, ModiaLang

makedocs(
  #modules  = [Modia],
  sitename = "ModiaLang",
  authors  = "Hilding Elmqvist (Mogram) and Martin Otter (DLR-SR)",
  format = Documenter.HTML(prettyurls = false),
  pages    = [
     "Home"      => "index.md",
  ]  
)
