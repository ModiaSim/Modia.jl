using Documenter, TinyModia, ModiaPlot

makedocs(
  #modules  = [TinyModia],
  sitename = "TinyModia",
  authors  = "Hilding Elmqvist (Mogram) and Martin Otter (DLR-SR)",
  format = Documenter.HTML(prettyurls = false),
  pages    = [
     "Home"      => "index.md",
	 "Tutorial"  => "Tutorial.md",
  	 "Functions" => "Functions.md",
     "Internal"  => "Internal.md"
  ]
)
