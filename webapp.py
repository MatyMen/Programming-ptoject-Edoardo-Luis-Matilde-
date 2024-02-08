from main import *
from flask import Flask, render_template, redirect, request, url_for
import pandas as pd
import os
import read_fasta




app = Flask(__name__)
app.myDNAClass = None
app.myRNAClass = None
app.myAAChainDict = None

@app.route("/", methods=['POST','GET'])
def home():


    return render_template("home.html")


@app.route('/setfile', methods=['POST', 'GET'])
def setfile():
    if 'file' not in request.files:
        return redirect(url_for('home'))
    file = request.files['file']
    # If the user does not select a file, the browser submits an
    # empty file without a filename.

    if file.filename == '':
        return redirect(url_for('home'))
    elif file.filename[-6:] != '.fasta' and file.filename[-3:] != '.fa' and file.filename[-4:] != '.fna':
        return redirect(url_for('home'))
    else:
        filename = "uploaded.fasta"
        os.makedirs(os.path.join(app.instance_path,"fastas"),exist_ok=True)
        os.makedirs(os.path.join( app.root_path,"static"), exist_ok=True) #this part of the code is used if static is lost, it makes possible to the software to work but css and explanation imeges are lost
        file.save(os.path.join(app.instance_path,"fastas", filename))
        myDataset=read_fasta.Dataset()
        myDataset.readfasta_set(os.path.join(app.instance_path, "fastas", 'uploaded.fasta'), True)

        app.myDNAClass = DNA(myDataset.get_data(),myDataset.get_DNA_RNA())

        return render_template('choice.html')

@app.route('/setfile1', methods=['POST', 'GET'])
def setfile1():
    if 'file' not in request.files:
        return redirect(url_for('home'))
    file = request.files['file']
    # If the user does not select a file, the browser submits an
    # empty file without a filename.

    if file.filename == '':
        return redirect(url_for('home'))
    elif file.filename[-6:] != '.fasta' and file.filename[-3:] != '.fa'  and file.filename[-4:] != '.fna' :
        return redirect(url_for('home'))
    else:
        filename = "uploaded.fasta"
        os.makedirs(os.path.join(app.instance_path,"fastas"),exist_ok=True)
        os.makedirs(os.path.join( app.root_path,"static"), exist_ok=True)
        file.save(os.path.join(app.instance_path,"fastas", filename))
        file_one = open((os.path.join(app.instance_path,"fastas", filename)), "a") #all fasta files that contains a protein miss the stop codon so the program add one
        file_one.write(".")
        file_one.close()
        myDataset = read_fasta.Dataset()
        myDataset.readfasta_set(os.path.join(app.instance_path, "fastas", 'uploaded.fasta'), False)

        app.myAAChainDict = DNA(myDataset.get_data(), myDataset.get_DNA_RNA())
        return redirect(url_for('translate_single'))



@app.route("/", methods=['POST', 'GET'])
def choice():
    return render_template("choice.html")






@app.route("/generate.html")
def generate():
    if app.myDNAClass == None:
        return redirect(url_for('home'))
    dna = app.myDNAClass.get_strand()
    dna_str = DNA.turn_in_str(dna)


    return render_template("generate.html", DNA=dna_str, NAME=app.myDNAClass.name[1:])

@app.route("/generate-.html")
def generateminus():
    if app.myDNAClass == None:
        return redirect(url_for('home'))
    if app.myDNAClass == None:
        return redirect(url_for('home'))

    dna = app.myDNAClass.produce_negative_strand()
    dna_str = DNA.turn_in_str(dna)






    return render_template("generate-.html", DNA=dna_str, NAME=app.myDNAClass.name[1:])


@app.route("/visualise.html")
def visualise():
    if app.myDNAClass == None:
        return redirect(url_for('home'))
    plt=Sequences.produce_graph(app.myDNAClass.strand.reset_index(drop=True))
    plt.savefig(os.path.join(app.root_path, "static",'graph.png'))
    data = DNA.gen_data(app.myDNAClass.get_strand())
    app.myRNAClass = app.myDNAClass.transcription()






    return render_template("visualise.html",DATA=data, NAME=app.myDNAClass.name[1:])


@app.route("/transcribe.html")
def transcribe():
    if app.myDNAClass == None:
        return redirect(url_for('home'))

    RNA_pos=app.myRNAClass.get_strand()
    RNA_neg=app.myRNAClass.produce_negative_strandRNA()
    RNA_pos_str=Sequences.turn_in_str(RNA_pos)
    RNA_neg_str=Sequences.turn_in_str(RNA_neg)




    return render_template("transcribe.html",RNA_POS=RNA_pos_str,RNA_NEG= RNA_neg_str,NAME=app.myDNAClass.name[1:] )

@app.route("/translate.html")
def translate():
    app.myAAChainDict = app.myRNAClass.translation()

    return render_template("translate.html")

@app.route("/translate1.html")
def translate1():
    if app.myDNAClass == None  :
        return redirect(url_for('home'))
    AAChain_strand=Sequences.turn_in_str(app.myAAChainDict[1].get_strand())

    return render_template("translate1.html",AACHAIN1= AAChain_strand,NAME=app.myDNAClass.name[1:])

@app.route("/translate2.html")
def translate2():
    if app.myDNAClass == None  :
        return redirect(url_for('home'))

    AAChain_strand=Sequences.turn_in_str(app.myAAChainDict[2].get_strand())

    return render_template("translate2.html",AACHAIN1= AAChain_strand,NAME=app.myDNAClass.name[1:])
@app.route("/translate3.html")
def translate3():
    if app.myDNAClass == None :
        return redirect(url_for('home'))

    AAChain_strand=Sequences.turn_in_str(app.myAAChainDict[3].get_strand())

    return render_template("translate3.html",AACHAIN1= AAChain_strand,NAME=app.myDNAClass.name[1:])
@app.route("/translate-1.html")
def translate_1():
    if app.myDNAClass == None  :
        return redirect(url_for('home'))

    AAChain_strand=Sequences.turn_in_str(app.myAAChainDict[-1].get_strand())

    return render_template("translate-1.html",AACHAIN1= AAChain_strand,NAME=app.myDNAClass.name[1:])
@app.route("/translate-2.html")
def translate_2():
    if app.myDNAClass == None  :
        return redirect(url_for('home'))

    AAChain_strand=Sequences.turn_in_str(app.myAAChainDict[-2].get_strand())

    return render_template("translate-2.html",AACHAIN1= AAChain_strand,NAME=app.myDNAClass.name[1:])
@app.route("/translate-3.html")
def translate_3():
    if app.myDNAClass == None  :
        return redirect(url_for('home'))

    AAChain_strand=Sequences.turn_in_str(app.myAAChainDict[-3].get_strand())

    return render_template("translate-3.html",AACHAIN1= AAChain_strand,NAME=app.myDNAClass.name[1:])
@app.route("/translate_single.html")
def translate_single():
    if app.myDNAClass == None and app.myAAChainDict==None :
        return redirect(url_for('home'))
    AAChain_strand=Sequences.turn_in_str(app.myAAChainDict.get_strand())
    plt = Sequences.produce_graph(app.myAAChainDict.get_strand())
    plt.savefig(os.path.join(app.root_path, "static", 'graphprot.png'))
    data = DNA.gen_data(app.myAAChainDict.get_strand())

    return render_template("translate_single.html",AACHAIN1= AAChain_strand,DATA=data,NAME=app.myAAChainDict.name[1:])

@app.route("/oligo.html")
def oligopeptides():
    if app.myDNAClass == None and app.myAAChainDict==None:
        return redirect(url_for('home'))
    dtorf_list = []
    for key in app.myAAChainDict.keys():
        dtorf_temp= app.myAAChainDict[key].get_oligos()
        dtorf_temp= dtorf_temp.assign(orf=key)
        dtorf_list.append(dtorf_temp)

    dtorffinal=pd.concat(dtorf_list)
    dtorffinal= dtorffinal.sort_values(by='length')
    Data=[]
    for index, row in dtorffinal.iterrows():
        curr_chain_str = ''
        for aa in row['aachains']:
            curr_chain_str = curr_chain_str + aa.get_letter()
        Data.append((curr_chain_str,row['orf'],index,row['length']))

    return render_template("oligo.html", DATA=Data,NAME=app.myDNAClass.name[1:])


@app.route("/proteins.html")
def proteins():
    if app.myDNAClass == None :
        return redirect(url_for('home'))
    dtorf_list = []
    for key in app.myAAChainDict.keys():
        dtorf_temp = app.myAAChainDict[key].get_proteins()

        dtorf_temp = dtorf_temp.assign(orf=key)
        dtorf_list.append(dtorf_temp)

    dtorffinalextreme = pd.concat(dtorf_list)
    dtorffinalextreme = dtorffinalextreme.sort_values(by='length')
    Data = []
    for index, row in dtorffinalextreme.iterrows():
        curr_chain_str = ''
        for aa in row['aachains']:
            curr_chain_str = curr_chain_str + aa.get_letter()
        Data.append((curr_chain_str, row['orf'], index, row['length']))

    return render_template("proteins.html", DATA=Data,NAME=app.myDNAClass.name[1:])



@app.route("/oligo12.html")
def singleoligopeptide():
    if app.myDNAClass == None :
        return redirect(url_for('home'))
    orf = request.args.get("orf")
    index = request.args.get("index")
    myAAChains = app.myAAChainDict[int(orf)]
    write1 = orf
    write2 = index
    write3 = myAAChains.get_single_prot_len(int(index))
    chain = myAAChains.get_single_aachain(int(index))
    plt = Sequences.produce_graph(chain)
    plt.savefig(os.path.join(app.root_path, "static", 'graphprot.png'))
    curr_chain_str=Sequences.turn_in_str(chain)
    data=Sequences.gen_data(chain)


    return render_template("oligo12.html",  WRITE1=write1, WRITE2=write2, WRITE3=write3, CHAIN=curr_chain_str,DATA=data)


@app.route("/protein12.html")
def singleprotein():
    if app.myDNAClass == None:
        return redirect(url_for('home'))
    orf = request.args.get("orf")
    index = request.args.get("index")
    myAAChains =app.myAAChainDict[int(orf)]
    write1 = orf
    write2 = index
    write3 = myAAChains.get_single_prot_len(int(index))
    chain = myAAChains.get_single_aachain(int(index))
    plt1 = Sequences.produce_graph(chain.reset_index(drop=True))
    plt1.savefig(os.path.join(app.root_path, "static", 'graphprot.png'))

    curr_chain_str = Sequences.turn_in_str(chain)
    data = Sequences.gen_data(chain)


    return render_template("protein12.html", WRITE1=write1, WRITE2=write2, WRITE3=write3, CHAIN=curr_chain_str,DATA=data )


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080, debug=False)
