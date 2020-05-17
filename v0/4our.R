library(shiny)
library(nglShiny)
library(htmlwidgets)
library(shinyBS)
library(yaml)
#----------------------------------------------------------------------------------------------------
tooltips <- yaml.load_file("tooltips.yaml")
for(i in 1:length(tooltips)) tooltips[[i]]$text <- paste(tooltips[[i]]$text, collapse=" ")
printf("length of tooltips read: %d", length(tooltips))
#----------------------------------------------------------------------------------------------------
printf <- function(...)print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#                 PAS 116–222 95 0.81
#                 GAF 234–432 185 0.73
#                 PHY 480–610 80 1.8
# PHY without hairpin 480–610 43 0.97
#             Hairpin 560–591
#----------------------------------------------------------------------------------------------------
# from https://www.uniprot.org/uniprot/P14713
components.4our=list(
    A=list(name="A",
           selection=":A",
           representation="cartoon",
           colorScheme="residueIndex",
           visible=FALSE),
    PAS=list(name="PAS",
                  selection="116-222 AND :A",
                  representation="cartoon",
                  colorScheme="residueIndex",
                  visible=TRUE),
    GAF=list(name="GAF",
                  selection="234-432 AND :A",
                  representation="cartoon",
                  colorScheme="residueIndex",
                  visible=FALSE),
    PHY=list(name="PHY",
                  selection="480-610 AND :A",
                  representation="cartoon",
                  colorScheme="residueIndex",
                  visible=TRUE),
    pyrrole.D=list(name="pyrrole.D",
                   selection="[2VO] AND :A and (.C01 OR .C02 OR .C03 OR .C04 OR .C05 OR .C06 OR .C42 OR .N41 OR .O43)",
                   representation="ball+stick",
                   colorScheme="element",
                   visible=TRUE),
    hairpinA=list(name="hairpinA",
                  selection="560-591 AND :A",
                  representation="cartoon",
                  colorScheme="residueIndex",
                  visible=FALSE),
    chromophoreA=list(name="chromophoreA",
                      selection="[2VO] AND :A",
                      representation="ball+stick",
                      colorScheme="element",
                      visible=FALSE)
    ) # components.4our



nglRepresentations = c('angle', 'axes', 'ball+stick', 'backbone', 'base', 'cartoon', 'contact',
                       'dihedral', 'distance', 'helixorient', 'licorice', 'hyperball', 'label',
                       'line', 'surface', 'point', 'ribbon', 'rocket', 'rope', 'spacefill', 'trace', 'unitcell',
                       'validation')
nglColorSchemes <- c('residueIndex', 'chainIndex', 'entityType', 'entityIndex')
defaultRepresentation <- "cartoon"
defaultColorScheme <- "residueIndex"
#----------------------------------------------------------------------------------------------------
# 1RQK, 3I4D: Photosynthetic reaction center from rhodobacter sphaeroides 2.4.1
# crambin, 1crn: https://bmcbiophys.biomedcentral.com/articles/10.1186/s13628-014-0008-0
addResourcePath("www", "www");

ui = shinyUI(fluidPage(

  tags$head(
    tags$style("#nglShiny_4our{height:90vh !important;}"),
    tags$link(rel="icon", href="data:;base64,iVBORw0KGgo=")
    ),

  includeCSS("phytochrome.css"),
  with(tooltips[[1]], bsTooltip(selector, text, location, options = list(container = "body", html=TRUE))),

  tabsetPanel(type = "tabs",
              tabPanel("Introduction",  includeHTML("intro.html")),
              tabPanel("PhyB",
                       sidebarLayout(
                           sidebarPanel(
                               actionButton("fitButton", "Fit"),
                                        #actionButton("defaultViewButton", "Defaults"),
                               actionButton("hideAllRepresentationsButton", "Hide All"),
                               actionButton("toggleAChainVisibilityButton", "Show All"),
                               br(), br(),
                                        #actionButton("toggleChromophoreBVisibilityButton", "Chromophore B"),
                                        #actionButton("toggleBChainVisibilityButton", ":B"),
                                        #actionButton("togglePASdomainVisibilityButton", "PAS"),
                                        #actionButton("toggleGAFdomainVisibilityButton", "GAF"),
                                        #br(); br()
                               h5("Domains"),
                               actionButton("togglePASVisibilityButton", "PAS"),
                               actionButton("toggleGAFVisibilityButton", "GAF"),
                               actionButton("togglePHYVisibilityButton", "PHY"),
                               actionButton("toggleHairpinAVisibilityButton", "PHY Hairpin"),

                               br(), br(),
                               h5("Prosthetic group"),
                               actionButton("toggleChromophoreAVisibilityButton", "Chromophore"),
                               actionButton("togglePyrroleDVisibilityButton", "Pyrrole D"),
                               width=2),
                           mainPanel(nglShinyOutput('nglShiny_4our'),width=10)
                           ),   # sidebarLayout
                       ), # PhyB tabPanel
              tabPanel("Notes", includeHTML("4our-notes.html")),
              tabPanel("Chromophore", includeHTML("chromophore.html")),
              tabPanel("Terms", includeHTML("terms.html")),
              tabPanel("Papers", includeHTML("papers.html"))
              ) # tabsetPanel
  ) # fluidPage
) # shinyUI
#----------------------------------------------------------------------------------------------------
server = function(input, output, session) {

   hideAll <- function(){
      components.4our$chromophoreA$visible <<- FALSE
      components.4our$chromophoreB$visible <<- FALSE
      components.4our$A$visible <<- FALSE
      components.4our$B$visible <<- FALSE
      components.4our$hairpinA$visible <<- FALSE
      components.4our$PAS$visible <<- FALSE
      components.4our$GAF$visible <<- FALSE
      components.4our$PHY$visible <<- FALSE
      components.4our$pyrrole.D$visible <<- FALSE

      setVisibility(session, "chromophoreA", FALSE)
      setVisibility(session, "chromophoreB", FALSE)
      setVisibility(session, "A", FALSE)
      setVisibility(session, "B", FALSE)
      setVisibility(session, "hairpinA", FALSE)
      setVisibility(session, "PAS", FALSE)
      setVisibility(session, "GAF", FALSE)
      setVisibility(session, "PHY", FALSE)
      setVisibility(session, "pyrrole.D", FALSE)

      }

  observeEvent(input$fitButton, ignoreInit=TRUE, {
     fit(session)
     })

  observeEvent(input$domainChooser, ignoreInit=TRUE, {
     chosenDomain = input$domainChooser
     printf("domains choice: %s", chosenDomain)
     residueRange <- switch(chosenDomain,
                            "helix001" = 7:19,
                            "helix002" = 23:30,
                            "sheet001" = 1:4,
                            "sheet002" = 32:35
                            )

     residue.string <- paste(residueRange, collapse=", ")
     printf("%s: %s", chosenDomain, residue.string)
     session$sendCustomMessage(type="select", message=list(residue.string))
       # 327 atoms, 46 residues
       # HELIX    1  H1 ILE A    7  PRO A   19  13/10 CONFORMATION RES 17,19       13
       # HELIX    2  H2 GLU A   23  THR A   30  1DISTORTED 3/10 AT RES 30           8
       # SHEET    1  S1 2 THR A   1  CYS A   4  0
       # SHEET    2  S1 2 CYS A  32  ILE A  35 -1
     })

  observeEvent(input$defaultViewButton, ignoreInit=TRUE, {
     session$sendCustomMessage(type="removeAllRepresentations", message=list())
     session$sendCustomMessage(type="setRepresentation", message=list(defaultRepresentation))
     session$sendCustomMessage(type="setColorScheme", message=list(defaultColorScheme))
     session$sendCustomMessage(type="fit", message=list())
     })


   observeEvent(input$showChromophoreAttachmentSiteButton, ignoreInit=TRUE, {
     repString <- "ball+stick"
     selectionString <- "24"
     session$sendCustomMessage(type="showSelection", message=list(representation=repString,
                                                                  selection=selectionString,
                                                                  name="chromophoreAttachment"))
     })

   observeEvent(input$toggleChromophoreAVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$chromophoreA$visible
     components.4our$chromophoreA$visible <<- newState
     components.4our$pyrrole.D$visible <<- newState
     setVisibility(session, "chromophoreA", newState)
     setVisibility(session, "pyrrole.D", newState)
     })

   observeEvent(input$toggleChromophoreBVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$chromophoreB$visible
     components.4our$chromophoreB$visible <<- newState
     setVisibility(session, "chromophoreB", newState)
     })


   observeEvent(input$toggleHairpinAVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$hairpinA$visible
     components.4our$hairpinA$visible <<- newState
     setVisibility(session, "hairpinA", newState)
     })

   observeEvent(input$togglePASVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$PAS$visible
     components.4our$PAS$visible <<- newState
     setVisibility(session, "PAS", newState)
     })

   observeEvent(input$toggleGAFVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$GAF$visible
     components.4our$GAF$visible <<- newState
     setVisibility(session, "GAF", newState)
     })

   observeEvent(input$togglePHYVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$PHY$visible
     components.4our$PHY$visible <<- newState
     setVisibility(session, "PHY", newState)
     })

   observeEvent(input$togglePyrroleDVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$pyrrole.D$visible
     components.4our$pyrrole.D$visible <<- newState
     setVisibility(session, "pyrrole.D", newState)
     })

   observeEvent(input$toggleAChainVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$A$visible
     components.4our$A$visible <<- newState
     setVisibility(session, "A", newState)
     })

   observeEvent(input$showAllButton, ignoreInit=TRUE, {
     hideAll()
     setVisibility(session, "A", newState)
     })

   observeEvent(input$toggleBChainVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$B$visible
     components.4our$B$visible <<- newState
     setVisibility(session, "B", newState)
     })

   observeEvent(input$toggleGAFdomainVisibilityButton, ignoreInit=TRUE, {
     newState <- !components.4our$gaf$visible
     components.4our$gaf$visible <<- newState
     setVisibility(session, "gaf", newState)
     })

   observeEvent(input$showCBDButton, ignoreInit=TRUE, {
     repString <- "cartoon"
     selectionString <- "1-321"
     colorScheme = "residueIndex"
     session$sendCustomMessage(type="showSelection", message=list(representation=repString,
                                                                  selection=selectionString,
                                                                  colorScheme=colorScheme,
                                                                  name="CBD"))
     })


   observeEvent(input$showCBD.PAS.Button, ignoreInit=TRUE, {
     repString <- "cartoon"
     selectionString <- "38-128"
     colorScheme = "residueIndex"
     session$sendCustomMessage(type="showSelection", message=list(representation=repString,
                                                                  selection=selectionString,
                                                                  colorScheme=colorScheme,
                                                                  name="PAS"))
     })

  observeEvent(input$showCBD.GAF.Button, ignoreInit=TRUE, {
     repString <- "cartoon"
     selectionString <- "129-321"
     colorScheme = "residueIndex"
     session$sendCustomMessage(type="showSelection", message=list(representation=repString,
                                                                  selection=selectionString,
                                                                  colorScheme=colorScheme,
                                                                  name="GAF"))
     })

  observeEvent(input$hideAllRepresentationsButton, ignoreInit=TRUE, {
     hideAll()
#      components.4our$chromophoreA$visible <<- FALSE
#      components.4our$chromophoreB$visible <<- FALSE
#      components.4our$A$visible <<- FALSE
#      components.4our$B$visible <<- FALSE
#      components.4our$hairpinA$visible <<- FALSE
#      components.4our$PAS$visible <<- FALSE
#      components.4our$GAF$visible <<- FALSE
#      components.4our$PHY$visible <<- FALSE
#      components.4our$pyrrole.D$visible <<- FALSE
#
#
#      setVisibility(session, "chromophoreA", FALSE)
#      setVisibility(session, "chromophoreB", FALSE)
#      setVisibility(session, "A", FALSE)
#      setVisibility(session, "B", FALSE)
#      setVisibility(session, "hairpinA", FALSE)
#      setVisibility(session, "PAS", FALSE)
#      setVisibility(session, "GAF", FALSE)
#      setVisibility(session, "PHY", FALSE)
#      setVisibility(session, "pyrrole.D", FALSE)
#      #updateSelectInput(session, "representationSelector", label=NULL, choices=NULL,  selected=defaultRepresentation)
#      #updateSelectInput(session, "colorSchemeSelector", label=NULL, choices=NULL,  selected=defaultColorScheme)
      })

  observeEvent(input$pdbSelector, ignoreInit=TRUE, {
     choice = input$pdbSelector
     printf("pdb: %s", choice)
     session$sendCustomMessage(type="setPDB", message=list(choice))
     updateSelectInput(session, "pdbSelector", label=NULL, choices=NULL,  selected=choice)
     })

  observeEvent(input$representationSelector, ignoreInit=TRUE, {
     choice = input$representationSelector;
     printf("rep: %s", choice)
     session$sendCustomMessage(type="setRepresentation", message=list(choice))
     updateSelectInput(session, "representationSelector", label=NULL, choices=NULL,  selected=choice)
     })

  observeEvent(input$colorSchemeSelector, ignoreInit=TRUE, {
     choice = input$colorSchemeSelector;
     printf("colorScheme: %s", choice)
     session$sendCustomMessage(type="setColorScheme", message=list(choice))
     updateSelectInput(session, "colorSchemeSelector", label=NULL, choices=NULL,  selected=choice)
     })

  output$value <- renderPrint({input$action})

  options.4our <- list(pdbID="4our", namedComponents=components.4our)

  output$nglShiny_4our <- renderNglShiny(
    nglShiny(options.4our, 300, 300)
    )

} # server
#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui=ui, server=server), port=5669)
#shinyApp(ui=ui, server=server)


