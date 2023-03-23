import gradio as gr
import pandas as pd
import pymysql
import requests
import io
import igraph
import zipfile

def get_protein_network(identifiers, species, required_score, network_flavor, network_type, hide_disconnected_nodes):
    # 蛋白列表
    identifiers = identifiers.split("\n")
    # 删除空值
    identifiers = list(filter(None, identifiers))
    # split 默认单引号改双引号
    identifiers = str(identifiers).replace("'", "\"")
    # 构造 html 字符串
    iframedoc = (
        f"""<!DOCTYPE html>
        <html>
          <head>
            <script type="text/javascript" src="https://string-db.org/javascript/combined_embedded_network_v2.0.4.js"></script>
            <script>
              function send_request() {{
                  getSTRING("https://string-db.org", {{
                      "species": "{species}",
                      "identifiers": {identifiers},
                      "required_score": "{required_score}",
                      "network_flavor": "{network_flavor}",
                      "network_type": "{network_type}",
                      "hide_disconnected_nodes": "{hide_disconnected_nodes}",
                      "caller_identity": "http://www.eyeseeworld.com/cbd/"
                  }})
              }}
            </script>
            <style>
              #stringEmbedded svg {{
                  display: block;
                  margin: auto;
              }}
            </style>
          </head>
          <body onload="send_request();">
            <p><b>Network:</b></p>
            <div id="stringEmbedded"></div>
          </body>
        </html>
        """
    )

    html = f"""
    <div class="iframe-container" style="position: relative;overflow: hidden;padding-top: 56.25%">
    <iframe id="myframe" style="position: absolute;width: 100%;height: 100%;top: 0;left: 0;margin:0 auto" name="result" allow="midi; geolocation; 
    microphone; camera; display-capture; encrypted-media;" sandbox="allow-modals allow-forms allow-scripts 
    allow-same-origin allow-popups allow-top-navigation-by-user-activation allow-downloads" allowfullscreen="" 
    allowpaymentrequest="" frameborder="0" srcdoc='{iframedoc}'></iframe></div>"""

    return html

def mysql_get_protein():
    # 连接 MySQL
    con = pymysql.connect(host='47.98.184.123', user='guest', password='guest_cbd', database='cbd_limina_top', charset='utf8mb4')
    # 定义 SQL 语句
    query = '''SELECT Biomarker, String_Name 
                  FROM biomarker 
                  WHERE Category="Protein" AND String_Name <>"" AND String_Name <>"NA" 
                  ORDER BY String_Name ASC'''
    # pandas 读取
    df = pd.read_sql(sql=query, con=con)
    value_list_tmp = df.values.tolist()
    value_list = ['   '.join(value[::-1]) for value in value_list_tmp]
    # 关闭连接
    con.close()
    return value_list

def choose_protein(proteins):
    symbol = [protein.split(" ")[0] for protein in proteins]
    return "\n".join(symbol)

def get_protein_from_file(file):
    protein_df = pd.read_csv(file.name, header=None, engine='python')
    symbol = protein_df[0].tolist()
    return "\n".join(symbol)

def clear_input():
    return None, None, None, None

def example_proteins(key):
    application = key.split()[2]
    # 连接 MySQL
    con = pymysql.connect(host='47.98.184.123', user='guest', password='guest_cbd', database='cbd_limina_top', charset='utf8mb4')
    # 定义 SQL 语句
    query = f'''SELECT String_Name 
                   FROM biomarker 
                   WHERE Category="Protein" AND (Application like "%{application}%") AND String_Name <>"" AND String_Name <>"NA" 
                   ORDER BY String_Name ASC'''
    # pandas 读取
    df = pd.read_sql(sql=query, con=con)
    # 去重
    value_list_tmp = set(df.values.flatten())
    # 转列表并按首字母排序
    value_list = list(value_list_tmp)
    value_list.sort()
    value_input = "\n".join(value_list)
    # 关闭连接
    con.close()
    # return protein_input.update(value=value_input)
    return value_input

def get_network_stats(identifiers, species):
    # 获取网络参数 tsv
    identifiers = identifiers.split("\n")
    identifiers = list(filter(None, identifiers))
    identifiers = "%0d".join(identifiers)
    url = f"https://string-db.org/api/tsv/ppi_enrichment?identifiers={identifiers}&species={species}"
    response = requests.get(url)
    data = response.text
    df = pd.read_csv(io.StringIO(data), sep="\t")
    # return df
    df_T = df.transpose().reset_index()
    df_T.columns = ["Parameters", "Value"]
    return df_T

def calculate_topo(identifiers, species, required_score, network_type):
    # 获取相互作用网络 json
    identifiers = identifiers.split("\n")
    identifiers = list(filter(None, identifiers))
    identifiers = "%0d".join(identifiers)
    url = f"https://string-db.org/api/json/network?identifiers={identifiers}&species={species}&required_score={required_score}&network_type={network_type}"
    response = requests.get(url)
    data = response.json()
    # 边元组的列表
    edgelist = []
    for i in data:
        edgelist.append((i["preferredName_A"], i["preferredName_B"]))
    edgelist = list(set(edgelist))

    # igraph 建立网络和计算
    g = igraph.Graph.TupleList(edgelist, directed = False)
    symbol = g.vs["name"]
    degree = g.degree(symbol)
    betweenness = g.betweenness(symbol)
    closeness = g.closeness(symbol)
    df = pd.DataFrame(list(zip(symbol, degree, betweenness, closeness)), 
                              columns = ["Symbol", "Degree", "Betweenness", "Closeness"])
    df = df.round({'Betweenness': 6, 'Closeness': 6})
    return df

def get_enrichment(identifiers, species):
    # 获取富集分析 tsv
    identifiers = identifiers.split("\n")
    identifiers = list(filter(None, identifiers))
    identifiers = "%0d".join(identifiers)
    url = f"https://string-db.org/api/tsv/enrichment?identifiers={identifiers}&species={species}"
    response = requests.get(url)
    data = response.text
    df = pd.read_csv(io.StringIO(data), sep="\t")

    # 处理不必要的列
    df["genes_in_background"] = df["number_of_genes"].astype(str) + '/' + df["number_of_genes_in_background"].astype(str)
    df.drop(["number_of_genes", "number_of_genes_in_background", "ncbiTaxonId", "preferredNames"], axis=1, inplace=True)
    # 重新指定列的顺序
    cols = list(df.columns)
    cols.remove("genes_in_background")
    cols.insert(2, "genes_in_background")
    df = df[cols]

    # 添加超链接
    url_dict = {
        "COMPARTMENTS": {
            "url_head": f"https://compartments.jensenlab.org/Entity?order=textmining,knowledge,predictions&knowledge=10&textmining=10&predictions=10&type1=-22&type2={species}&id1=GO:",
            "color": "#ffb142"
        },
        "TISSUES": {
            "url_head": f"https://tissues.jensenlab.org/Entity?order=textmining,knowledge,experiments&knowledge=10&experiments=10&textmining=10&type1=-25&type2={species}&id1=",
            "color": "#B33771"
        },
        "DISEASES": {
            "url_head": f"https://diseases.jensenlab.org/Entity?order=textmining,knowledge,experiments&textmining=10&knowledge=10&experiments=10&type1=-26&type2={species}&id1=",
            "color": "#eb2f06"
        },
        "Process": {
            "url_head": f"http://amigo.geneontology.org/amigo/term/",
            "color": "#33d9b2"
        },
        "Component": {
            "url_head": f"http://amigo.geneontology.org/amigo/term/",
            "color": "#ff5252"
        },
        "Function": {
            "url_head": f"http://amigo.geneontology.org/amigo/term/",
            "color": "#34ace0"
        },
        "Keyword": {
            "url_head": f"https://www.uniprot.org/keywords/",
            "color": "black"
        },
        "KEGG": {
            "url_head": f"https://www.kegg.jp/kegg-bin/show_pathway?",
            "color": "#f78fb3"
        },
        "SMART": {
            "url_head": f"http://smart.embl-heidelberg.de/smart/do_annotation.pl?DOMAIN=",
            "color": "#BDC581"
        },
        "InterPro": {
            "url_head": f"https://www.ebi.ac.uk/interpro/entry/",
            "color": "#2daec1"
        },
        "Pfam": {
            "url_head": f"https://pfam.xfam.org/family/",
            "color": "#8854d0"
        },
        "PMID": {
            "url_head": f"https://pubmed.ncbi.nlm.nih.gov/",
            "color": "#20558a"
        },
        "RCTM": {
            "url_head": f"https://reactome.org/content/detail/R-",
            "color": "#778beb"
        },
        "WikiPathways": {
            "url_head": f"https://www.wikipathways.org/index.php/Pathway:",
            "color": "#a5b1c2"
        },
        "HPO": {
            "url_head": f"https://monarchinitiative.org/phenotype/",
            "color": "#4b6584"
        },
        "NetworkNeighborAL": {
            "url_head": f"https://string-db.org/cgi/network?input_query_species={species}&network_term_id=",
            "color": "#e17055"
        }
    }

    line = -1
    for term in df["term"]:
        line += 1
        category = df.iloc[line, 0]
        url_head = url_dict[category]["url_head"]
        color = url_dict[category]["color"]
        if category in ("COMPARTMENTS", "PMID"):
            term_text = term.split(':')[1]
        else:
            term_text = term
        df.iloc[line, 1] = f'<a href="{url_head}{term_text}" target="_blank" style="color: {color};font-weight:bold;">{term}</a>'

    return df

def visible(identifiers):
    identifiers = identifiers.split("\n")
    identifiers = list(filter(None, identifiers))
    if len(identifiers) > 0:
        return topo_para.update(visible=True), network_stats.update(visible=True), enrichment.update(visible=True), button_download.update(visible=True)

def download_all(identifiers, species, required_score, network_type):
    # 网络图像
    url1 = f"https://string-db.org/api/image/network?identifiers={identifiers}&species={species}"
    # 相互作用网络文件
    url2 = f"https://string-db.org/api/tsv/network?identifiers={identifiers}&species={species}&required_score={required_score}&network_type={network_type}"
    # 富集分析
    url3 = f"https://string-db.org/api/tsv/enrichment?identifiers={identifiers}&species={species}"

    # url 列表
    urls = [url1, url2, url3]
    # 定义文件名
    filenames = ["network.png", "network_interactions.tsv", "network_enrichment.tsv"]

    # 创建一个空的zip文件对象
    z = zipfile.ZipFile("cbd_string_analysis.zip", "w")
    for url, filename in zip(urls, filenames):
        r = requests.get(url)
        # 创建一个类似于文件的对象，用于存储内容
        file_obj = io.BytesIO(r.content)
        # 将内容写入zip文件
        z.writestr(filename, file_obj.read())
    # 关闭zip文件对象
    z.close()
    return "cbd_string_analysis.zip"

def download_visible(identifiers, species, required_score, network_type):
    file = download_all(identifiers, species, required_score, network_type)
    if file != None:
        return file_download.update(visible=True)

example = f"""
    <div class="iframe-container" style="position: relative;overflow: hidden;padding-top: 56.25%">
    <iframe id="myframe" style="position: absolute;width: 100%;height: 100%;top: 0;left: 0;margin:0 auto" name="result" allow="midi; geolocation; 
    microphone; camera; display-capture; encrypted-media;" sandbox="allow-modals allow-forms allow-scripts 
    allow-same-origin allow-popups allow-top-navigation-by-user-activation allow-downloads" allowfullscreen="" 
    allowpaymentrequest="" frameborder="0" 
    srcdoc='<!DOCTYPE html>
                <html>
                  <head>
                    <script type="text/javascript" src="https://string-db.org/javascript/combined_embedded_network_v2.0.4.js"></script>
                    <script>
                      function send_request() {{
                          getSTRING("https://string-db.org", {{
                              "species": "9606",
                              "identifiers": ["TP53"],
                              "required_score": "400",
                              "network_flavor": "confidence",
                              "network_type": "functional",
                              "hide_disconnected_nodes": "0"
                          }})
                      }}
                    </script>
                    <style>
                      #stringEmbedded svg {{
                          display: block;
                          margin: auto;
                      }}
                    </style>
                  </head>
                  <body onload="send_request();">
                    <h3 style="text-align: justify;">
                    The Explore page is designed to allow the user to select or enter any biomarker (or other proteins) of interest 
                    to obtain a network of interactions between them and to allow the user to query the network for relevant information 
                    such as topological parameters, enrichment information, etc. The available parameters are already listed on the left side.</h3>
                    <h3 style="text-align: justify;">
                    In addition, if you enter only one protein, you will obtain its interaction with 10 other proteins (based on STRING scores). 
                    Here is an example for input protein "TP53":</h3>
                    <div id="stringEmbedded"></div>
                  </body>
                </html>'>
    </iframe>
    </div>"""

css = """
:root {
    --radius-lg: 5px;
    --color-background-primary: #f9fbff;
    --color-background-tertiary: #f9fbff;
}

body {
    background-color: #f9fbff;
}

/*滚动条*/
*::-webkit-scrollbar {
    width: 7px;
    height: 7px;
}

*::-webkit-scrollbar-thumb {
    background: #014371;
}

/* 兼容Firefox、IE*/
* {
    scrollbar-width: 10px;
    scrollbar-base-color: green;
    scrollbar-track-color: red;
    scrollbar-arrow-color: blue;
}

#checkbox {
    height: 270px;
    overflow: auto !important;
}

[data-testid="checkbox-group"] {
    display: block !important;
}

#submit-button {
    background: #00d4ad;
    color: white;
    border: 2px solid #00d4ad;
}

#download-button {
    background: #00add4;
    color: white;
    border: 2px solid #00add4;
}

#clear-button:active,
#submit-button:active,
#download-button:active {
    box-shadow: inset 0 5px 10px #e1e5ee;
}

.table-wrap.svelte-1sckpv1.svelte-1sckpv1.svelte-1sckpv1 {
    border: none;
}

#network-stats table, 
#topo-para table, 
#enrichment table {
    border-collapse: collapse;
    box-shadow: 3px 0px 5px #ddd;
    background-color: white;
    text-align: center;
    overflow: hidden;
    word-break: break-word;
    margin: auto;
}

#topo-para thead tr, 
#enrichment thead tr {
    width: calc(100% - 7px);
}

#topo-para tbody {
    display: block;
    overflow-y: scroll;
    max-height: 300px;
    width: 100%;
}

#network-stats thead, 
#topo-para thead, 
#enrichment thead {
    box-shadow: 0 5px 10px #e1e5ee;
}

#network-stats th, 
#topo-para th, 
#enrichment th {
    padding: 0.7rem 1rem;
    text-transform: uppercase;
    letter-spacing: 0.1rem;
    font-size: 0.9rem;
    font-weight: 900;
    text-align: center;
}

#network-stats td, 
#topo-para td, 
#enrichment td {
    padding: 0.5rem 0.7rem;
    box-shadow: 3px 0px 5px #ddd;
}

#enrichment tbody {
    display: block;
    overflow-y: scroll;
    max-height: 800px;
    width: 100%;
}

#network-stats tbody tr:nth-child(even), 
#topo-para tbody tr:nth-child(even), 
#enrichment tbody tr:nth-child(even) {
    background: #d3e4f5;
}

#topo-para tr, 
#enrichment tr {
    box-sizing: border-box;
    table-layout: fixed;
    display: table;
    width: 100%;
}

#enrichment tbody tr:hover {
    background: #E0EEE0;
    transform: translateY(-3px);
    transition: all .3s;
}

#enrichment tbody a {
    position: relative;
    text-decoration: none;
}

#enrichment tbody a::before {
    content: '';
    position: absolute;
    width: 0;
    height: 2px;
    bottom: 0;
    left: 0;
    background-color: #216a9d;
    transition: width 0.2s ease-in-out;
}

#enrichment tbody a:hover::before {
    width: 100%;
}

#network-stats th span,
#network-stats td span,
#topo-para th span,
#topo-para td span {
    overflow: hidden;
    word-wrap: break-word;
    white-space: normal;
}

#enrichment th span,
#enrichment td span {
    min-width: 100px;
    overflow: hidden;
    word-wrap: break-word;
    word-break: break-all;
    white-space: normal;
}
"""

# 使用gr.Blocks()创建和组合组件
with gr.Blocks(css=css) as demo:
    with gr.Row():
        with gr.Column(scale=1):
            # 创建输入组件
            with gr.Tab("Proteins in CBD2"):
                protein_choose = gr.CheckboxGroup(choices=mysql_get_protein(), label='Proteins', info='Proteins in CBD2', elem_id='checkbox')
            with gr.Tab("Proteins from file"):
                protein_file = gr.File(file_types=['text', '.csv'])
            protein_input = gr.TextArea(max_lines = 10, label='Protein Name', placeholder='Three ways to submit your query for protein(s):\n1. Type here;\n2. Select the protein(s) in CBD2 above;\n3. Upload your protein file(.txt or .csv) above.\nNote: One protein per line.')
            species_input = gr.Dropdown(choices=['9606', '10090', '10116'], value='9606', label='Species')
            slider_input = gr.Slider(minimum=0, maximum=1000, value=400, label='Score Threshold')
            network_flavor_input = gr.Dropdown(choices=['confidence', 'evidence', 'actions'], value='confidence', label='Network Flavor')
            network_type_input = gr.Dropdown(choices=['functional', 'physical'], value='functional', label='Network Type')
            hide_disconnected_input = gr.Dropdown(choices=['0', '1'], value='0', label='Hide Disconnected Nodes')
            
            # 创建按钮组件
            with gr.Row():
                button_clear = gr.Button("Clear", elem_id="clear-button")
                button_input = gr.Button("Submit", elem_id="submit-button")

            # 临时组件
            cache = gr.Textbox(visible=False)

            # 下崽按钮
            # button_download = gr.Button("Download Analysis", visible=False, elem_id="download-button")

            # 网络参数输出列表
            network_stats = gr.Dataframe(label="Network Stats:", elem_id="network-stats", visible=False)
            # topo_para = gr.Dataframe(label="Topological Parameters:", elem_id="topo-para", visible=False)
            # gr.update(visible=True)
            
        with gr.Column(scale=3.25, min_width=985):
            # HTML 输出组件
            html_output = gr.HTML(value=example)
            with gr.Row():
                # 例子组件
                gr.Examples(examples=[
                                                  "Biomarkers for Diagnosis in CBD2", 
                                                  "Biomarkers for Prognosis in CBD2", 
                                                  "Biomarkers for Treatment in CBD2"
                                                  ], 
                                  label="Try Examples", 
                                  inputs=cache, 
                                  outputs=[protein_input], 
                                  fn=example_proteins, 
                                  cache_examples=True, 
                                  run_on_click=True)
                # 下崽组件
                button_download = gr.Button("Download Analysis", visible=False, elem_id="download-button").style(full_width=False)
            file_download = gr.File(label="Download Analysis", visible=False, elem_id="download-file")
            # 网络参数组件
            # network_stats = gr.Dataframe(label="Network Stats:", elem_id="network-stats", visible=False)
            topo_para = gr.Dataframe(label="Topological Parameters:", elem_id="topo-para", visible=False)
            enrichment = gr.Dataframe(label="Functional Enrichments:", 
                                                 datatype=["str", "html", "str", "str", "number", "number", "str"], 
                                                 elem_id="enrichment", visible=False)

    # 按钮监听
    # 清空输入
    button_clear.click(fn=clear_input, inputs=None, outputs=[protein_choose, protein_file, protein_input, cache], queue=False)
    # 获得网络 svg
    button_input.click(fn=get_protein_network, 
                             inputs=[protein_input, species_input, slider_input, network_flavor_input, network_type_input, hide_disconnected_input], 
                             outputs=html_output)
    # 获得网络参数
    button_input.click(fn=get_network_stats, 
                             inputs=[protein_input, species_input], 
                             outputs=network_stats)
    # 计算拓扑参数
    button_input.click(fn=calculate_topo, 
                             inputs=[protein_input, species_input, slider_input, network_type_input], 
                             outputs=topo_para)
    # 获得功能富集结果
    button_input.click(fn=get_enrichment, 
                             inputs=[protein_input, species_input], 
                             outputs=enrichment)
    # 显示列表
    button_input.click(fn=visible, inputs=protein_input, outputs=[topo_para, network_stats, enrichment, button_download])
    # 下载分析结果
    button_download.click(fn=download_all, inputs=[protein_input, species_input, slider_input, network_type_input], outputs=file_download)
    button_download.click(fn=download_visible, inputs=[protein_input, species_input, slider_input, network_type_input], outputs=file_download)

    # 获取蛋白写入文本框
    protein_choose.change(fn=choose_protein, inputs=protein_choose, outputs=protein_input)
    protein_file.change(fn=get_protein_from_file, inputs=protein_file, outputs=protein_input)

# 在本地运行demo对象
demo.launch(share=False)
