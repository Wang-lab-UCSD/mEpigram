"""usage: julia motifscanner_B.jl p-valuecutoff memefile positivefile DSnegativefile output"""
function score(kmer,PWM,background)
    #calculate log-odds scores for the matches
    pseudo=0.001
    s=0
    for i =1:endof(kmer)
        #s+=log(PWM[kmer[i]][i]+pseudo)
        s+=log(PWM[kmer[i]][i]+pseudo) - log(background[kmer[i]]+pseudo)
    end
    return s
end
function load_motifs(filename)
    file=open(filename)
    seq=split(strip(readstring(file)),"MOTIF")
    seq=seq[2:end];
    motifs=Dict()
    #ans=Dict{AbstractString,Any}[]
    for s=1:endof(seq)
        t=split(strip(seq[s]),"\n")
        motifs[t[1]]=t[3:end]
    end
    for m in keys(motifs)
        #println(m)
        tdict=Dict('A'=>Float64[],'C'=>Float64[],'G'=>Float64[],'T'=>Float64[],'E'=>Float64[])
        for pos=1:endof(motifs[m])
            #println(pos)
            tmp=split(strip(motifs[m][pos]),"\t")
            #print(tmp)
            push!(tdict['A'],float(tmp[1]))
            push!(tdict['C'],float(tmp[2]))
            push!(tdict['G'],float(tmp[3]))
            push!(tdict['T'],float(tmp[4]))
            push!(tdict['E'],float(tmp[5]))
        end 
        motifs[m]=tdict
        #ans[m]=tdict
    
    end
    return motifs
end 

function revcompl(word,transdict)
    #return word
    revComp=Char[]
    index=1
    while index <=endof(word)
        if word[index] == 'E' #if is methylated C
            #println("is E",index)
            if index < endof(word) && word[index+1]=='G' #if G is next to it
                append!(revComp,['G','E'])
                #println("found")
                index+=2   
            else
                #println("adding G, condition not met")
                append!(revComp,['G'])
                index+=1
            end 
        else
            #println("adding",transdict[word[index]])
            append!(revComp,[transdict[word[index]]])
            index+=1  
        end 
    end 
    return revComp[end:-1:1]
end

function scoreCutoff(pwm,DsSeqs,pvaluecutoff,bg)
    """
    Calculate the scores for all of the kmers in the shuffled sequences, get the distribution then decide the score cutoff
    corresponding to the p-value cutoff
    """
    k=length(pwm['A'])
    scores=Float64[]
    for seq in DsSeqs
        if length(seq)<2*k  || 'N' in seq
            continue
        else
            tmp=collect(seq)
            for i=1:length(tmp)-k+1
                kmer=tmp[i:i+k-1]
                revcomlkmer=revcompl(kmer,transdict)
                push!(scores,max(score(kmer,pwm,bg),score(revcomlkmer,pwm,bg))) #get the maximum score between the revcompl and the forward
            end
        end
    end 
    num=length(scores)
    scores=sort(scores,rev=true)
    scorecutoff=scores[Int(floor(num*pvaluecutoff))]*0.999
    #println( "score cutoff is",scorecutoff)
    return scorecutoff,scores[1],scores[end]
end


function scoreCutoff_new(pwm,DsSeqs,pvaluecutoff,bg) #dont use this one
    """
    Calculate the scores for all of the kmers in the shuffled sequences, get the distribution then decide the score cutoff
    corresponding to the p-value cutoff
    """
    k=length(pwm['A'])
    scores=Float64[]
    totalscore=0
    for seq in DsSeqs
        if length(seq)<2*k  || 'N' in seq
            continue
        else
            tmp=collect(seq)
            for i=1:length(tmp)-k+1
                kmer=tmp[i:i+k-1]
                revcomlkmer=revcompl(kmer,transdict)
                s=max(score(kmer,pwm,bg),score(revcomlkmer,pwm,bg))
                totalscore+=s
                push!(scores,s) #get the maximum score between the revcompl and the forward
            end
        end
    end 
    num=length(scores)
    mean=totalscore/num
    std=var(scores)^0.5
    println(string(num)*" scores, mean is "*string(mean)*", std is "*string(std))
    inname="tmp."*string(mean)*".in.txt"
    outname="tmp."*string(mean)*".out.txt"
    target=open(inname,"w")
    write(target,string(mean)*"\t"*string(std)*"\t"*string(pvaluecutoff))
    close(target)
    
    run(`python ./calcscorecutoff.py $inname $outname`)
    scorecutoff=float(strip(readstring(open(outname))))
    return scorecutoff,scores[1],scores[end]
end
function var(array)
    mean=sum(array)/length(array)
    variance=0
    for number in array
        variance+=(number-mean)^2
    end 
    return variance/length(array)
end 
            
function main()
	#bg=Dict('A'=> 0.205, 'T' => 0.205, 'C' => 0.2874, 'G' => 0.295, 'E' => 0.0076)
    p_cutoff = float(ARGS[1])
    motiffile = ARGS[2]
    positiveregions = ARGS[3]
    negativeregions = ARGS[4]
    output = ARGS[5]

    bgfile =ARGS[6]
    seq=split(strip(readstring(open(bgfile))),'\n')
    bg=Dict()
    for line in seq
        println(line)
        tmp=split(strip(line),'\t')
        bg[collect((tmp[1]))[1]]=float(tmp[2])
    end 
    #println(bg)
    
    #load the background
    println("Scanning "*positiveregions*" with "*motiffile*" at pvalue "*string(p_cutoff))
    motifs=load_motifs(motiffile)
    #println(typeof(motifs))
    trainnegset = split(strip(readstring(open(negativeregions))),'>')[2:end]
    posset = split(strip(readstring(open(positiveregions))),'>')[2:end]
    
    global transdict=Dict('A'=>'T','C'=>'G','G'=>'C','T'=>'A')
    trainnegseqs = AbstractString[] 
    #trainnegseqs_names=[]
    for seq in trainnegset
        try
            tmp=split(strip(seq),"\n")
            push!(trainnegseqs,tmp[2])
        catch e
            println(e)
            continue
        end 
    end 
    
    posseqs=AbstractString[]
    posseqs_names=AbstractString[]
    for seq in posset
        try
            tmp=split(strip(seq),"\n")
            if length(tmp)!=2
            	continue
            end
            push!(posseqs,tmp[2])
            push!(posseqs_names,tmp[1])
        catch e
            println(e)
            continue
        end 
    end 
    motifnames=keys(motifs)
    #motifnames=["98_4.821_0-6-6-12-23-23-29-29-29-23-23-17-6-6-0-0_8_29"]
    #print(typeof(motifs))
    target=open(output,"w")
    for m in motifnames
        pwm=motifs[m]
        println(m)
        write(target,"MOTIF\t"*m*"\n")
        trainnegscorecutoff,trainnegmaxscore,trainnegminscore=scoreCutoff(pwm,trainnegseqs,p_cutoff,bg)
        println("max score ",trainnegmaxscore)
        println("min score ",trainnegminscore)
        println("score cutoff ",trainnegscorecutoff)
        
        #scanning positive seqs
        counts=0
        k=length(pwm['A'])
        for i=1:length(posseqs)
            seq = posseqs[i]
            seqname = posseqs_names[i]
            if length(seq)<2*k  || 'N' in seq
                continue
            else 
                tmp = collect(seq)
                for i=1:length(tmp)-k
                    kmer = tmp[i:i+k-1]
                    revcomlkmer=revcompl(kmer,transdict)
                    ts1=score(kmer,pwm,bg)
                    ts2=score(revcomlkmer,pwm,bg)
                    s=max(ts1,ts2)
                    if s>trainnegscorecutoff
                        counts+=1
                        if ts1 >=ts2
                        	write(target,seqname*"\t"*string(i)*"\t"*join(kmer)*"\t"*"+"*"\t"*string(s)*"\n")
                        else
                        	write(target,seqname*"\t"*string(i)*"\t"*join(kmer)*"\t"*"-"*"\t"*string(s)*"\n")
                        end
                    end 
                end
            end 
        end
        println("TOTAL\t"*string(counts)*" matches\n")
        write(target,"TOTAL\t"*string(counts)*"\n")
    end 
    close(target)
end 

@time(main())