module ModiaLogging

logTranslationDefault = false # if logging of translation should be made (in a file in the folder testresults or in the terminal window).
global logOnFile # if logging should made on file or in terminal window
global logName = "Modia_log"

# ----------------------------------------------

export openLogModia, logModia, loglnModia, closeLogModia, logFileModia, setDefaultLogName, resetTestStatus, setTestStatus, increaseLogCategory, printTestStatus

global logTranslation = logTranslationDefault

@static if VERSION < v"0.7.0-DEV.2005"
    global defaultOutput = STDOUT
else
    global defaultOutput = stdout
end

global logFileModia = defaultOutput
global defaultLogName
setDefaultLogName(name) = global defaultLogName = name

function setOptions(options) 
    global logTranslation = logTranslationDefault
    if haskey(options, :logTranslation)
        global logTranslation = options[:logTranslation]
        @show logTranslation
        delete!(options, :logTranslation)
    end

    global logOnFile = true
    if haskey(options, :logOnFile)
        global logOnFile = options[:logOnFile]
        @show logOnFile
        delete!(options, :logOnFile)
    end
  
    global logName = defaultLogName
    if haskey(options, :logName)
        global logName = options[:logName]
        @show logName
        delete!(options, :logName)
    end
end

function openLogModia() 
    if logOnFile && logTranslation
        logDirectory = homedir() * "/" * "ModiaResults"

        if !isdir(logDirectory)
            mkdir(homedir() * "/" * "ModiaResults")
        end

        logFileName = homedir() * "/" * "ModiaResults/" * logName * ".txt"
        global logFileModia = open(logFileName, "w")
        println("Log file: ", logFileName)
    else
        global logFileModia = defaultOutput
    end
    # write(logFileModia, "\r\n") # Trying to enable viewing log in Notepad
end

logModia(args...) = if logTranslation; print(logFileModia, args...) else end
loglnModia(args...) = if logTranslation; println(logFileModia, args...) else end
closeLogModia() = if logOnFile && logTranslation; close(logFileModia) else end

global nOK = 0
global nNOTOK = 0
global logCategories = Dict()

function resetTestStatus()
    global nOK = 0
    global nNOTOK = 0
    global logCategories = Dict()
    return
end

function setTestStatus(OK)
    if OK
        global nOK += 1
    else
        global nNOTOK += 1
    end
end

function increaseLogCategory(category)
    if haskey(logCategories, category)
        global logCategories[category] += 1
    else
        global logCategories[category] = 1
    end
end

@static if VERSION < v"0.7.0-DEV.2005"
    printstyled(s, i; bold=false, color=:black) = print_with_color(color, s, i, bold=bold)
end  

function printTestStatus()
    global nOK
    global nNOTOK
    global logCategories
    println()
    println("\n----------------------\n")
    println()
    printstyled("Number of simulations OK    : ", nOK, bold=true, color=:green); println()
    printstyled("Number of simulations NOT OK: ", nNOTOK, bold=true, color=:red); println()
    println()
    println("Log category statistics:")
    for (cat, count) in logCategories
        println(cat, ": ", count)
    end    
    println("\n----------------------\n")
    println()
    resetTestStatus()
end

end
