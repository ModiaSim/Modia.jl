using Documenter, Modia, Modia.SignalTables, SignalTablesInterface_PyPlot

makedocs(
  #modules  = [Modia],
  sitename = "Modia",
  authors  = "Hilding Elmqvist (Mogram) and Martin Otter (DLR-SR)",
  format = Documenter.HTML(prettyurls = false),
  pages    = [
     "Home"      => "index.md",
	 "Tutorial"  => [
       "tutorial/Tutorial.md"
       "tutorial/GettingStarted.md"
       "tutorial/Modeling.md"
       "tutorial/Simulation.md"
       "tutorial/FloatingPointTypes.md"
       "tutorial/Appendix.md"
      ],
  	 "Functions" => "Functions.md",
     "Internal"  => "Internal.md"
  ]  
)
