#!/usr/bin/env python3
#
# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0
#

# Python import(s).
import argparse
import csv

# ROOT import(s).
import ROOT

def runOnInputs(inputs, func):
    '''Function running another function on every row in all input CSVs

    Argument(s):
       inputs  -- List of CSV file name, platform name pairs
       func    -- The function-like-object to run on every CSV row
    '''

    # Loop over the files.
    for input in inputs:
        # Open this particular file.
        with open(input[0], 'r') as csvFile:
            csvReader = csv.DictReader(csvFile)
            # Process each row of the file.
            for row in csvReader:
                # Execute the received function on this row.
                func(row)
                pass
            pass
        pass
    return

def getSampleNames(inputs):
    '''Get the union of "sample names" from the inputs

    Argument(s):
       inputs  -- List of CSV file name, platform name pairs

    Return:
       A list of unique sample names
    '''

    # The result list.
    result = []

    def collectDirectory(row):
        '''Collect the directory names from the CSV rows

        Argument(s):
           row -- A CSV row coming from a CSV DictReader
        '''
        # Use the "result" variable from the parent scope.
        nonlocal result
        # If this directory was not seen yet, add it to the result.
        dirname = row['directory']
        if not dirname in result:
            result.append(row['directory'])
            pass
        return

    # Collect the directory names.
    runOnInputs(inputs, collectDirectory)

    # Return the list.
    return result

def getMaxThreads(inputs):
    '''Get the maximum number of CPU threads for which (a) measurement(s) exist

    Argument(s):
       inputs  -- List of CSV file name, platform name pairs

    Return:
       The maximum number of threads found in the input files
    '''

    # The result value.
    result = 1

    def updateThread(row):
        '''Find the maximum number of threads in the CSV rows

        Argument(s):
           row -- A CSV row coming from a CSV DictReader
        '''
        # Use the "result" variable from the parent scope.
        nonlocal result
        # If we see a value larger than what we saw before, let's use it.
        threads = int(row['threads'])
        if threads > result:
            result = threads
            pass
        return

    # Find the maximum number of threads.
    runOnInputs(inputs, updateThread)

    # Return the value.
    return result

def getMinMaxThroughput(inputs, sample = ''):
    '''Get the minimum/maximum throughput values

    Argument(s):
       inputs  -- List of CSV file name, platform name pairs
       sample  -- A specific sample to look at, or all samples if left empty

    Return:
       A pair of minimum-maximum throughput values
    '''

    # The result pair.
    result = [1E10, 0.]

    def findMinMax(row):
        '''Find the minimum and maximum throughput value in the input files

        Argument(s):
           row -- A CSV row coming from a CSV DictReader
        '''
        # Use the variable(s) from the parent scope.
        nonlocal sample
        nonlocal result
        # Skip the current row if it's not coming from the appropriate
        # directory/sample.
        if sample != '' and row['directory'] != sample:
            return
        # Update the min/max values.
        throughput = getThroughput(row)
        if throughput < result[0]:
            result[0] = throughput
            pass
        if throughput > result[1]:
            result[1] = throughput
            pass
        return

    # Find the minimum-maximum throughput values.
    runOnInputs(inputs, findMinMax)

    # Return the minimum-maximum pair.
    return result

def getThroughput(csvRow):
    '''Get the throughput in Hz for a given CSV row

    Argument(s):
       csvRow  -- A CSV row coming from a CSV DictReader

    Returns:
       The throughput value in Hz
    '''

    return (float(csvRow['processed_events']) /
            (float(csvRow['processing_time']) * 1E-9))

def configureProfile(profile, platform, graphicsOptions):
    '''Configure a TProfile object with the appropriate properties

    Argument(s):
       profile          -- The TProfile object to configure
       platform         -- The platform that the profile represents
       graphicsOptions  -- The graphics options to use for the profile

    Return:
       The configured TProfile object
    '''

    # Set all (relevant) properties on the profile object.
    profile.SetMarkerStyle(graphicsOptions[platform]['MarkerStyle'])
    profile.SetMarkerColor(graphicsOptions[platform]['MarkerColor'])
    profile.SetMarkerSize (graphicsOptions[platform]['MarkerSize'])
    profile.SetLineColor(graphicsOptions[platform]['LineColor'])

    # Return the configured TProfile.
    return profile

def drawPerMuSamplePlots(inputs, samples, canvas, output, maxThreads,
                         graphicsOptions):
    '''Function drawing the per-mu-sample plots

    Argument(s):
       inputs          -- List of CSV file name, platform name pairs
       samples         -- Names of the (mu) samples to make plots for
       canvas          -- The TCanvas to draw on
       output          -- Output PDF file name
       maxThreads      -- Maximum number of threads to make the plots for
       graphicsOptions -- The graphics options to use for the profile
    '''

    # Dummy counter used in giving unique names to the profile objects.
    globalCounter = 1

    # Flag helping with the PDF writing.
    firstPage = True

    # The result list with maximum througputs per sample, per platform.
    result = []

    # Loop over the samples.
    for sample in samples:

        # Helper variable for the draw option(s).
        drawOptions = 'E1 P'
        # Helper counter for the inputs.
        inputCounter = 0
        # List of temporary profile objects. These need to stay alive until the
        # canvas is printed.
        profiles = []

        # Get minimum-maximum throghput values for this sample.
        min_max_throughput = getMinMaxThroughput(inputs, sample)

        # The legend for all of the platforms used.
        legend = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)

        # Loop over the inputs.
        for input in inputs:

            # Create a TProfile for the results.
            profile = configureProfile(ROOT.TProfile(
                'throughput_%i' % globalCounter,
                'Throughput for "%s";CPU Threads;Events/Second' % sample,
                maxThreads, 0.5, maxThreads + 0.5), input[1], graphicsOptions)
            profile.GetYaxis().SetRangeUser(0., min_max_throughput[1] * 1.2)
            profiles.append(profile)

            # Add the profile to the legend.
            legend.AddEntry(profile, input[1], 'P')

            # Increment the counter(s).
            globalCounter += 1
            inputCounter += 1

            # Get the maximum throughput for this sample and platform.
            max_throughput = 0.

            # Open the input file.
            with open(input[0], 'r') as csvFile:
                csvReader = csv.DictReader(csvFile)
                # Process each row of the file.
                for row in csvReader:
                    # Skip measurements from samples other than the current one.
                    if row['directory'] != sample:
                        continue
                    # Add the result to the plot.
                    throughput = getThroughput(row)
                    profile.Fill(float(row['threads']), throughput)
                    if throughput > max_throughput:
                        max_throughput = throughput
                        pass
                    pass
                pass

            # Save the max throughput.
            result.append([sample, input[1], max_throughput])

            # Draw the profile on the canvas.
            profile.Draw(drawOptions)
            drawOptions = 'SAME %s' % drawOptions
            pass

        # Draw the legend.
        legend.Draw()

        # Save the page with all the profiles.
        if firstPage:
            canvas.Print('%s(' % output, 'pdf')
            firstPage = False
        else:
            canvas.Print(output, 'pdf')
            pass
        pass

    # Return the max throughputs for each sample and platform.
    return result

def drawMaxThroughputPlot(inputs, max_throughputs, samples, canvas, output,
                          graphicsOptions):
    '''Draw a single plot with the maximum throughputs on each platform

    Argument(s):
       inputs          -- List of CSV file name, platform name pairs
       max_throughputs -- The throughputs collected in drawPerMuSamplePlots()
       samples         -- Names of the (mu) samples to make plots for
       canvas          -- The TCanvas to draw on
       output          -- Output PDF file name
       graphicsOptions -- The graphics options to use for the profile
    '''

    # Collect the platform names, and the minimum-maximum values.
    platforms = []
    for thr in max_throughputs:
        if not thr[1] in platforms:
            platforms.append(thr[1])
            pass
        pass

    # Get the absolute minimum-maximum throughput values.
    min_max_throughput = getMinMaxThroughput(inputs)

    # The legend for all of the platforms used.
    legend = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)

    # List of temporary profile objects. These need to stay alive until the
    # canvas is printed.
    profiles = []

    # Helper counter.
    counter = 0

    # Helper variable for the draw option(s).
    drawOptions = 'E1 P L'

    # Display the maximum throughputs in logarithmic scale.
    canvas.SetLogy()

    # Create one profile per platform.
    for platform in platforms:

        # Create the profile.
        profile = configureProfile(ROOT.TProfile(
            'max_throughput_profile_%i' % counter,
            'Maximum Throughput;;Events/Second',
            len(samples), 0.5, len(samples) + 0.5), platform, graphicsOptions)
        for i in range(len(samples)):
            profile.GetXaxis().SetBinLabel(i + 1, samples[i])
            pass
        profile.GetYaxis().SetRangeUser(min_max_throughput[0] * 0.2,
                                        min_max_throughput[1] * 5.0)
        profiles.append(profile)
        counter += 1

        # Add the profile to the legend.
        legend.AddEntry(profile, platform, 'P')

        # Fill it with data from the previous maximum collection.
        for thr in max_throughputs:

            # Skip other platforms.
            if thr[1] != platform:
                continue

            # Fill the profile.
            profile.Fill(thr[0], thr[2])
            pass

        # Draw the profile on the canvas.
        profile.Draw(drawOptions)
        drawOptions = 'SAME %s' % drawOptions
        pass

    # Draw the legend.
    legend.Draw()

    # Save the page with all the profiles.
    canvas.Print('%s)' % output, 'pdf')
    return

def main():
    '''C(++)-style main function
    '''

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(
        description='Multi-Threaded Throghput Plotter')
    parser.add_argument('csv_file', nargs='+',
                        help='CSV file name:Platform name')
    parser.add_argument('-o', '--output', help='Output PDF file name',
                        default='output.pdf')
    parser.add_argument('-t', '--max-threads',
                        help='Maximum number of CPU threads to plot for',
                        type=int, dest='max_threads', default=-1)
    args = parser.parse_args()

    # Parse the CSV file argument(s).
    inputs = []
    for arg in args.csv_file:
        inputs.append(arg.split(':', 1))
        if len(inputs[-1]) != 2:
            raise RuntimeError(
                'Could not parse "%s" as "<filename>:<platform>"' % arg)
        pass

    # Get the union of "sample names" from all the inputs.
    samples = getSampleNames(inputs)

    # Figure out how many threads to show the results for.
    maxThreads = args.max_threads
    if maxThreads <= 0:
        maxThreads = getMaxThreads(inputs)
        pass

    # Associate ROOT graphical properties with all of the platforms.
    graphicsOptions = {}
    counter = 0
    for input in inputs:
        # Choose a color.
        color = counter + 1
        # Skip yellow...
        if color == 5:
            color += 1
            pass
        # Save the options.
        graphicsOptions[input[1]] = {'MarkerColor': color,
                                    'MarkerStyle': counter + 20,
                                     'MarkerSize' : 1,
                                     'LineColor'  : color}
        counter += 1
        pass

    # Create a canvas to draw on.
    canvas = ROOT.TCanvas('canvas', 'Throughput Canvas', 800, 600)
    canvas.cd()

    # Draw the per-mu-sample througput plots. The return value is a list
    # of best througput values per mu value for each platform.
    maxThroughputs = \
        drawPerMuSamplePlots(inputs, samples, canvas, args.output, maxThreads,
                             graphicsOptions)

    # Draw the maximum throughput values for each mu value, for each platform.
    drawMaxThroughputPlot(inputs, maxThroughputs, samples, canvas, args.output,
                          graphicsOptions)

    # Return gracefully.
    return 0

if __name__ == '__main__':
    # Set ROOT into batch mode.
    ROOT.gROOT.SetBatch()
    # Set some global ROOT style option(s).
    ROOT.gStyle.SetOptStat(False)
    # Execute the main() function.
    import sys
    sys.exit(main())
