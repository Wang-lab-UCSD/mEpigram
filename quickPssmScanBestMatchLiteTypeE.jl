"""usage: julia quickPssmScanBestMatchLiteTypeE.jl memefile file.faa tagname resuldir bgComposition"""
function score(kmer,PWM,background)
    #calculate log-odds scores for the matches
    pseudo=0.001
    s=0
    for i =1:endof(kmer)
        char = kmer[i]
        s+=log(PWM[char][i]+pseudo) - log(background[char]+pseudo)
    end
    return s
end

function convertPWMtoMatrixE(PWM)
    """ convert the typeE PWM dict into a 2D matrix"""
    alphabet = ['A','C','G','T','E']
    alsize = length(alphabet)
    k = length(PWM['A'])
    matrix = Array{Float64}(alsize, k)    
    for pos=1:k
        for i=1:alsize
            char = alphabet[i]
            matrix[i,pos] = PWM[char][pos]
        end
    end
    return matrix
end

function score_w_matrix(kmer,matrix,bgmatrix)
    #calculate log-odds scores for the matches
    pseudo=0.001
    s=0
    for i=1:endof(kmer)
        #println(i)
        char = parse(Int,(kmer[i]))
        #println(char)
        s+=log(matrix[char,i]+pseudo) - log(bgmatrix[char]+pseudo)
    end
    return s
end

function seqtoNum(seq,transdict)
    """Convert DNA sequence to seuqnce of numbers"""
    newstr = "" 
    for i=1:endof(seq)
        newstr*= transdict[seq[i]]
    end
    return newstr 
end 


function load_motifs(filename)
    file=open(filename)
    seq=split(strip(readstring(file)),"MOTIF")
    seq=seq[2:end];
    motifs=Dict()
    #ans=Dict{ASCIIString,Any}[]
    for s=1:endof(seq)
        t=split(strip(seq[s]),"\n")
        motifs[t[1]]=t[3:end] #skip the first 2 lines of the meme 
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
    return join(revComp[end:-1:1])
end


function main()
    """This function takes in a memefile, a set of sequences, and output the best score for each sequence by each motif in each separate file"""
    motiffile = ARGS[1]
    sequencefile = ARGS[2]
    tag = ARGS[3]
    resultdir = ARGS[4] 
    bgfile =ARGS[5]

    seq=split(strip(readstring(open(bgfile))),'\n')
    bg=Dict()
    for line in seq
        tmp=split(strip(line),'\t')
        bg[collect((tmp[1]))[1]]=float(tmp[2])
    end
    alphabet = ['A','C','G','T','E'] #type E
    bgmatrix = Array{Float64}(length(alphabet))
    for i=1:length(alphabet)
        bgmatrix[i] = bg[alphabet[i]]
    end

    motifs=load_motifs(motiffile)
    #println("Loaded motifs")
    
    transdict=Dict('A'=>'T','C'=>'G','G'=>'C','T'=>'A', 'N'=>'N')
    transdict2=Dict('A'=>"1",'C'=>"2",'G'=>"3",'T'=>"4",'E'=>"5", 'F'=>"6",'N'=>"7")

    motifnames=keys(motifs)
    for m in motifnames
        pwm=motifs[m]
        pwmmatrix = convertPWMtoMatrixE(pwm)
        println("Scanning "*m)
        target=open(resultdir*"/"*m*"."*tag*".scanned.txt","w")
        write(target,"MOTIF\t"*m*"\n")
        counts=0
        k = length(pwm['A'])
        f = open(sequencefile)
        header =""
        for line in eachline(f)
            if line[1]=='>'
                header=strip(line)
            else
                seq=uppercase(strip(line))
                seqname=header
                revcomplseq=revcompl(seq,transdict)
                sequencestoscan=[seq;revcomplseq] #add them all together in a same list
                bestscore=-100
                for sequence in sequencestoscan
                    if length(sequence)<k
                        continue
                    else
                        sequence = seqtoNum(sequence, transdict2)#convert to num
                        for i=1:length(sequence)-k
                            kmer = sequence[i:i+k-1]
                            if '7' in kmer #kmer is a list of digits, 7 means  'N'
                                continue
                            else
                                #s=score(kmer,pwm,bg)
                                s = score_w_matrix(kmer,pwmmatrix,bgmatrix)
                                if s > bestscore
                                    bestscore=s
                                end
                            end
                        end
                    end 
                end
                write(target,seqname*"\t"*string(bestscore)*"\n")
            end
        end
        close(target)
    end 
end

@time(main())