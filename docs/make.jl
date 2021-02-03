using Documenter, TinyModia

makedocs(
  #modules  = [TinyModia],
  sitename = "TinyModia",
  authors  = "Hilding Elmqvist (Mogram) and Martin Otter (DLR-SR)",
  format = Documenter.HTML(prettyurls = false),
  pages    = [
     "Home"      => "index.md",
	 "TinyModia" => "TinyModia.md",
  	 "Functions" => "Functions.md",
     "Internal"  => "Internal.md"
  ]
)
