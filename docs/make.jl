using Documenter, Modia, ModiaLang, ModiaPlot

makedocs(
  #modules  = [Modia],
  sitename = "Modia",
  authors  = "Hilding Elmqvist (Mogram) and Martin Otter (DLR-SR)",
  format = Documenter.HTML(prettyurls = false),
  pages    = [
     "Home"      => "index.md",
	 "Modia Tutorial"  => Any[
                              "tutorial/GettingStarted.md",
                              "tutorial/Modeling.md",
                              "tutorial/Simulation.md",
                              "tutorial/FloatingPointTypes.md",
                              "tutorial/Appendix.md"
                          ],
  	 "Functions" => "Functions.md",
     "Internal"  => "Internal.md"
  ]
)
