from django.shortcuts import render
from django.http import HttpResponse
from django.template import loader
import pandas as pd
import geopandas as gpd
import pyproj
import os
from plotly.offline import plot
import plotly.express as px
import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
from dna_features_viewer import BiopythonTranslator
from matplotlib import use
from matplotlib import pyplot as plt
import json
from plotly.subplots import make_subplots

use('agg')

color_dict = {
    "Adenovirus": "#65E0C8",
    "Bocavirus": "#D2D2E3",
    "Influenza A": "#6875D9",
    "Influenza B": "#E158D3",
    "Metapneumovirus": "royalblue",
    "Parainfluenza 2": "beige",
    "Parainfluenza 3": "#848877",
    "Parainfluenza 4": "#FDBF6F",
    "Polyomavirus": "#D9664D",
    "RSV - A/B": "thistle",
    "SARS-CoV-2" : "#E9D588",
    "Rhino- / Enterovirus": "#74457B",
    "coronavirus 229E": "pink",
    "coronavirus HKU1": "#BAE6B0",
    "coronavirus NL63": "#A968E1",
    "coronavirus OC43": "darkblue",
    "Mixed": "darkorange",
    "Negative": "lightgrey"
}

def pileup_plot(coverage_file, annotation_file, out_file, figsize=(20, 5), height_ratios=[4, 1.5]):
    coverage=pd.read_csv(coverage_file)
    fig, (ax1, ax2) = plt.subplots(
        2, 1,  figsize=figsize, sharex=True, gridspec_kw={"height_ratios": height_ratios}
    )
    ax1.plot(coverage['idx'], coverage['mean'])
    ax1.fill_between(coverage['idx'], coverage['ci_lower'], coverage['ci_upper'], color='b', alpha=.15)
    ax1.set_ylabel("Sequencing Depth", fontsize = 10)
    ax1.set_title("Sequencing Depth" ,fontsize = 15, loc='left' ,pad=20)
    ax1.set_yscale("log")
    ax1.set_ylim(ymin=1)
    ax1.get_xaxis().set_visible(False)
    #ax1.axhline(y = 5, color = 'b', linestyle = '--', label = "DP5")
    ax1.axhline(y = 10, color = 'b', linestyle = '-', label = "DP10")
    #ax1.axhline(y = 20, color = 'b', linestyle = ':', label = "DP20")
    ax1.legend(bbox_to_anchor = (1, 1), loc = 'upper left')
    
    graphic_record = BiopythonTranslator().translate_record(annotation_file)
    graphic_record.plot(ax=ax2, strand_in_label_threshold=4, with_ruler=True)
    ax2.set_xlabel("Position")
    ax2.get_yaxis().set_visible(False)
    
    fig.savefig(out_file, bbox_inches='tight')

def home(request):
    
    df = pd.read_csv('dashboard/static/example_data.csv')

    df['ent_date']= pd.to_datetime(df['ent_date'])
    df['week']= df['ent_date'].dt.strftime('%Y/%W')


    df['strain_name'][df['pseudonymized_id'].duplicated(keep=False)] = 'Mixed'
    df['strains_PCR'] = df['strains_PCR'].fillna('Negative')
    df['strains_PCR'] = [s if len(s.split(', '))<=1 else 'Mixed' for s in df['strains_PCR'].astype(str)]

    df = df.drop_duplicates(subset='pseudonymized_id')
    
    no_samples = len(df)
    '''
    grouped_df = df.groupby(['week','strain_name']).agg({"pseudonymized_id": ["count"], "match_PCR_sequencing": ["sum"]}).reset_index()
    grouped_df.columns = ['week','strain_name','count', 'match']
    grouped_df['match'] = grouped_df['match'].astype(str) + '/' + grouped_df['count'].astype(str)
    #size().reset_index(name='count')
    grouped_df['week'] = grouped_df['week'] + '/1'
    grouped_df['week']= pd.to_datetime(grouped_df['week'], format='%Y/%W/%w')

    fig = px.bar(grouped_df, x="week", y="count", color="strain_name",
            hover_data=['count', 'match'], barmode = 'stack',
            labels={
                     "week": "Week",
                     "count": "No. infections",
                     "strain_name": "Strain",
                     "match": "Match PCR"
                 },
                color_discrete_map= color_dict,
                category_orders = {"strain_name": color_dict.keys()},
                template = "simple_white",
                title='Sequencing'
            )
    #fig.update_xaxes(categoryorder='array', categoryarray= grouped_df['week'])
    fig.update_yaxes(showgrid=True)
    fig.update_xaxes(tickformat = '%Y/%W')

    #bar_seq = plot(fig, output_type="div")
    fig.write_html('dashboard/static/histogram_seq.html')

    grouped_df = df.groupby(['week','strains_PCR']).agg({"pseudonymized_id": ["count"], "match_PCR_sequencing": ["sum"]}).reset_index()
    grouped_df.columns = ['week','strain_PCR','count', 'match']
    grouped_df['match'] = grouped_df['match'].astype(str) + '/' + grouped_df['count'].astype(str)
    grouped_df['week'] = grouped_df['week'] + '/1'
    grouped_df['week']= pd.to_datetime(grouped_df['week'], format='%Y/%W/%w')

    fig = px.bar(grouped_df, x="week", y="count", color="strain_PCR",
            hover_data=['count', "match"], barmode = 'stack',
            labels={
                     "week": "Week",
                     "count": "No. infections",
                     "strain_PCR": "Strain",
                     "match": "Match sequencing"
                 },
            color_discrete_map= color_dict,
            category_orders = {"strain_PCR": color_dict.keys()},
            template = "simple_white",
            title = 'PCR'
            )
    fig.update_xaxes(tickformat = '%Y/%W')
    fig.update_yaxes(showgrid=True)

    #bar_pcr = plot(fig, output_type="div")
    fig.write_html('dashboard/static/histogram_pcr.html')
    
    f="dashboard/static/swissBOUNDARIES3D_1_5_TLM_KANTONSGEBIET.shx"
    shapes = gpd.read_file(f, engine="pyogrio")

    shapes = shapes.to_crs('EPSG:4326')

    shapes['canton'] = ['GV', 'TG' ,'VS', 'AG', 'SZ', 'ZH' ,'OW', 'FR', 'GL', 'UR', 'NW', 'SO', 'AR', 'JU', 'GR', 'VD', 'LU', 'TI', 'ZG', 'BL', 'SG', 'SH', 'BE', 'BS', 'NE', 'AI']
    shapes=shapes[['geometry', 'canton', 'NAME']]

    shapes["geometry"] = (
        shapes.to_crs(shapes.estimate_utm_crs()).simplify(5000).to_crs(shapes.crs)
    )

    shapes.to_file('dashboard/static/swiss_cantons.geojson', driver='GeoJSON')
    
    
    geojson = gpd.read_file('dashboard/static/swiss_cantons.geojson')


    weeks = df['week'].unique()
    weeks.sort()

    cantons = np.array([[i]*len(weeks) for i in geojson.canton]).flatten()
    weeks = list(weeks) * len(geojson)

    data = pd.DataFrame({'canton': cantons, 'week': weeks})


    no_infections = df.groupby(['prescriber', 'week']).size().reset_index(name='count')
    no_infections.columns = ['canton', 'week', 'No. infections']
    
    data = data.merge(no_infections, how='left', on=['canton', 'week'])
    data['No. infections'] = data['No. infections'].fillna(0)

    data = geojson.merge(data, how='right', on= 'canton')
   
    fig = px.choropleth_mapbox(data,
                       geojson=data.geometry,
                       locations=data.index,
                       hover_name = 'NAME',
                       zoom = 5.5,
                       center = {"lat": 46.8, "lon": 8.5},
                       color='No. infections',
                       color_continuous_scale='viridis',
                       animation_frame='week',
                       mapbox_style = "carto-positron"
                       )
    texttrace = go.Scattermapbox(
        lat=geojson.geometry.centroid.y,
        lon=geojson.geometry.centroid.x,
        text=geojson["canton"],
        textfont={"color":"white","size":10},
        mode="text",
        name="canton",
        hoverinfo='skip'
    )

    fig.add_trace(texttrace)

    fig.update_layout(
    coloraxis_colorbar=dict(
        tickvals=list(range(int(max(data['No. infections']))+1)),
    )
)


    fig.write_html('dashboard/static/map.html')

    #map_plt = plot(fig, output_type="div")
    '''
    
    return render(request, 'home.html', {'no_samples': no_samples }) # 'barplot_seq': bar_seq, 'barplot_pcr': bar_pcr, 

def strain_view(request, strain_name):

    df = pd.read_csv('dashboard/static/example_data.csv')

    df['ent_date']= pd.to_datetime(df['ent_date'])
    df['week'] = df['ent_date'].map(lambda dt: dt.strftime('%Y/%W'))

    match strain_name:
        case 'Enterovirus':
            strain_name_long = 'Rhino- / Enterovirus'
        case 'RSV':
            strain_name_long = 'RSV - A/B'
        case _:
            strain_name_long= strain_name.replace('_' ,' ')

    pcr_positive = df.apply(lambda row: strain_name_long in str(row['strains_PCR']).split(', '), axis=1)
    
    df_pcr = df[pcr_positive].drop_duplicates(subset='pseudonymized_id')
    df = df[df['strain_name']==strain_name_long]
    no_samples = len(df)
    substrains = df['substrain_name'].unique().tolist()


    if not os.path.isfile('dashboard/static/barplots/'+ strain_name + '_seq.html'):

        if(strain_name not in ['Bocavirus', 'Parainfluenza_4', 'Polyomavirus']):

            grouped_df = df_pcr.groupby('week').agg({"pseudonymized_id": ["count"], "match_PCR_sequencing": ["sum"]}).reset_index()
            grouped_df.columns = ['week', 'count', 'match']
            grouped_df['match'] = grouped_df['match'].astype(str) + '/' + grouped_df['count'].astype(str)
            grouped_df['week'] = grouped_df['week'] + '/1'
            grouped_df['week']= pd.to_datetime(grouped_df['week'], format='%Y/%W/%w')

            fig = px.bar(grouped_df, x="week", y="count", 
                    labels={
                             "week": "Week",
                             "count": "No. infections",
                             "match": "Match sequencing"
                         },
                    hover_data=['count', "match"],
                    template = "simple_white")
            
            fig.update_xaxes(tickformat = '%Y/%W')
            fig.update_yaxes(showgrid=True)
            fig.update_layout(legend=dict(xanchor = "left",  x = 0, yanchor = "bottom", y = 1), yaxis=dict(tickvals=list(range(max(grouped_df['count'])+3)))) 

            fig.write_html('dashboard/static/barplots/'+ strain_name + '_pcr.html')


        if len(substrains) > 1:
            '''
            if strain_name == 'Enterovirus':
                df['substrain_name'] = list(map({"Human rhinovirus A89": "Rhinovirus A", "Human enterovirus C109 isolate NICA08-4327": "Enterovirus C"}.get, df['substrain_name']))
            else:
                df['substrain_name'] = list(map({"Respiratory syncytial virus (type A)": "RSV A", "Human Respiratory syncytial virus 9320 (type B)": "RSV B"}.get, df['substrain_name']))
            '''

            grouped_df = df.groupby(['week', 'substrain_name']).agg({"pseudonymized_id": ["count"], "match_PCR_sequencing": ["sum"]}).reset_index()
            grouped_df.columns = ['week', 'substrain_name', 'count', 'match']
            grouped_df['match'] = grouped_df['match'].astype(str) + '/' + grouped_df['count'].astype(str)
            grouped_df['week'] = grouped_df['week'] + '/1'
            grouped_df['week']= pd.to_datetime(grouped_df['week'], format='%Y/%W/%w')

            fig = px.bar(grouped_df, x="week", y="count", color="substrain_name",
                barmode = 'stack',
                labels={
                         "week": "Week",
                         "count": "No. infections",
                         "substrain_name": "Substrain",
                         "match": "Match PCR"
                        },
                hover_data=['substrain_name', 'count', 'match'],
                template = "simple_white"
                )
        
        else:

            grouped_df = df.groupby('week').agg({"pseudonymized_id": ["count"], "match_PCR_sequencing": ["sum"]}).reset_index()
            grouped_df.columns = ['week', 'count', 'match']
            grouped_df['match'] = grouped_df['match'].astype(str) + '/' + grouped_df['count'].astype(str)
            grouped_df['week'] = grouped_df['week'] + '/1'
            grouped_df['week']= pd.to_datetime(grouped_df['week'], format='%Y/%W/%w')

            fig = px.bar(grouped_df, x="week", y="count", 
                    labels={
                             "week": "Week",
                             "count": "No. infections",
                             "match": "Match PCR"
                         },
                    hover_data=['count', "match"],
                    template = "simple_white")
        
        fig.update_xaxes(tickformat = '%Y/%W')
        fig.update_yaxes(showgrid=True)
        fig.update_layout(legend=dict(xanchor = "left",  x = 0, yanchor = "bottom", y = 1), yaxis=dict(tickvals=list(range(max(grouped_df['count'])+3)))) 

        fig.write_html('dashboard/static/barplots/'+ strain_name + '_seq.html')

        #bar_plt = plot(fig, output_type="div")

    geojson = gpd.read_file('dashboard/static/swiss_cantons.geojson')


    weeks = df['week'].unique()
    weeks.sort()

    cantons = np.array([[i]*len(weeks) for i in geojson.canton]).flatten()
    weeks = list(weeks) * len(geojson)

    data = pd.DataFrame({'canton': cantons, 'week': weeks})


    no_infections = df.groupby(['prescriber', 'week']).size().reset_index(name='count')
    no_infections.columns = ['canton', 'week', 'No. infections']
    
    data = data.merge(no_infections, how='left', on=['canton', 'week'])
    data['No. infections'] = data['No. infections'].fillna(0)

    data = geojson.merge(data, how='right', on= 'canton')

    if len(substrains) > 1:
        substrain_df = pd.pivot_table(df, index=['prescriber','week'], columns='substrain_name', values='pseudonymized_id', aggfunc='count').reset_index()
        substrain_df.rename(columns={ 'prescriber': "canton" }, inplace = True)
        substrains = df['substrain_name'].unique()
        data = data.merge(substrain_df, how='left', on=['canton', 'week'])
        data = data.fillna(0)

        fig = px.choropleth_mapbox(data,
                       geojson=data.geometry,
                       locations=data.index,
                       hover_name = 'NAME',
                       hover_data = substrains,
                       zoom = 5.5,
                       center = {"lat": 46.8, "lon": 8.5},
                       color='No. infections',
                       color_continuous_scale='viridis',
                       animation_frame='week',
                       mapbox_style = "carto-positron"
                       )

    else:
        fig = px.choropleth_mapbox(data,
                           geojson=data.geometry,
                           locations=data.index,
                           hover_name = 'NAME',
                           zoom = 5.5,
                           center = {"lat": 46.8, "lon": 8.5},
                           color='No. infections',
                           color_continuous_scale='viridis',
                           animation_frame='week',
                           mapbox_style = "carto-positron"
                           )
    texttrace = go.Scattermapbox(
        lat=geojson.geometry.centroid.y,
        lon=geojson.geometry.centroid.x,
        text=geojson["canton"],
        textfont={"color":"white","size":10},
        mode="text",
        name="canton",
        hoverinfo='skip'
    )

    fig.add_trace(texttrace)
    fig.update_layout(
        coloraxis_colorbar=dict(
            tickvals=list(range(int(max(data['No. infections']))+1)),
        )
    )
    map_plt = plot(fig, output_type="div")

    match strain_name_long:
        case 'Influenza A':
            if not os.path.isfile('dashboard/static/pileup/Influenza_A_PB1_all.png'):
                
                segments = ['PB2','PB1','PA','HA','NP', 'NA','MP','NS']

                for seg in segments:
                    for substrain in ['all', 'H1N1', 'H3N2']:
                        pileup_plot('dashboard/static/pileup/Influenza_A_'+seg+'_'+substrain+'.csv', 'dashboard/static/annotations/Influenza_A_'+seg+'.gb', 'dashboard/static/pileup/Influenza_A_'+seg+'_'+substrain+'.png')
                        
        case 'Influenza B':
            if not os.path.isfile('dashboard/static/pileup/Influenza_B_PB1.png'):
                
                segments = ['PB2','PB1','HA','NP', 'NA','MP','NS']

                for seg in segments:
                    pileup_plot('dashboard/static/pileup/Influenza_B_'+seg+'_all.csv', 'dashboard/static/annotations/Influenza_B_'+seg+'.gb', 'dashboard/static/pileup/Influenza_B_'+seg+'_all.png')
                    
        case _:
            if strain_name == 'Adenovirus': 
                figsize = (20, 6.5)
                height_ratios = [4, 2.5]
            else:
                figsize = (20, 5)
                height_ratios = [4, 1]

            pileup_plot('dashboard/static/pileup/'+strain_name+'_all.csv', 'dashboard/static/annotations/'+strain_name+'.gb', 'dashboard/static/pileup/'+strain_name+'_all.png', figsize=figsize, height_ratios=height_ratios)

            if len(substrains) > 1:
                for substrain in substrains:
                    pileup_plot('dashboard/static/pileup/'+strain_name+'_'+substrain.replace(' ' ,'_')+'.csv', 'dashboard/static/annotations/'+strain_name+'.gb', 'dashboard/static/pileup/'+strain_name+'_'+substrain.replace(' ' ,'_')+'.png', figsize=figsize, height_ratios=height_ratios)
            elif no_samples <= 25:
                coverage=pd.read_csv('dashboard/static/pileup/'+strain_name + '_all_indiv.csv', index_col=0)
                fig, (ax1, ax2) = plt.subplots(
                    2, 1,  figsize=figsize, sharex=True, gridspec_kw={"height_ratios": height_ratios}
                )
                coverage.plot(ax=ax1)
                ax1.set_ylabel("Sequencing Depth", fontsize = 10)
                ax1.set_title("Sequencing Depth" ,fontsize = 15, loc='left' ,pad=20)
                ax1.set_yscale("log")
                ax1.set_ylim(ymin=1)
                ax1.get_xaxis().set_visible(False)
                #ax1.axhline(y = 5, color = 'b', linestyle = '--', label = "DP5")
                ax1.axhline(y = 10, color = 'b', linestyle = '-', label = "DP10")
                #ax1.axhline(y = 20, color = 'b', linestyle = ':', label = "DP20")
                ax1.legend(title='Sample', bbox_to_anchor=(1, 1), loc='upper left')
                                
                graphic_record = BiopythonTranslator().translate_record('dashboard/static/annotations/'+strain_name+'.gb')
                graphic_record.plot(ax=ax2, strand_in_label_threshold=4, with_ruler=True)
                ax2.set_xlabel("Position")
                ax2.get_yaxis().set_visible(False)
                
                fig.savefig('dashboard/static/pileup/'+strain_name+'_all_indiv.png', bbox_inches='tight')


    '''
    fig = go.Figure([

        go.Scatter(
                name='Avg. depth',
                x=coverage['idx'],
                y=coverage['mean'],
                mode='lines',
                line=dict(color='rgb(31, 119, 180)')
            ),
        go.Scatter(
                x=list(coverage['idx'])+list(coverage['idx'][::-1]), 
                y=list(coverage['ci_upper'])+list(coverage['ci_lower'][::-1]),
                fill='toself',
                fillcolor='rgba(0,100,80,0.2)',
                line=dict(color='rgba(255,255,255,0)'),
                hoverinfo='skip',
                showlegend=False,
                name='95% CI'
            )
    ])

    fig.update_layout(
        xaxis_title='Position',
        yaxis_title='Sequencing depth',
        title='Sequencing depth',
        hovermode='x',
        yaxis_range=[0,np.log10(max(coverage['ci_upper']))]
    )
    fig.update_yaxes(type="log")

    pileup_plt = plot(fig, output_type="div")
    '''

    return render(request, "strain.html", {"strain": strain_name, "strain_name": strain_name_long, 'map': map_plt , 'no_samples': no_samples, 'substrains': substrains }) #'barplot': bar_plt,

def mixed_view(request):

    df = pd.read_csv('dashboard/static/example_data.csv')

    df = df[df['pseudonymized_id'].duplicated(keep=False)][['pseudonymized_id', 'strain_name', 'substrain_name']]

    mat = pd.pivot_table(df, index='pseudonymized_id', columns='strain_name', values='substrain_name', aggfunc='count', fill_value=0).astype(int)
    co_infections = mat.T.dot(mat)

    co_infections_np = co_infections.values.astype(float)
    co_infections_np[np.triu_indices(len(co_infections_np))] = np.nan

    no_co_infections = int(np.nansum(co_infections_np))

    fig = go.Figure(go.Heatmap(z=co_infections_np[1:,:-1][::-1], x=co_infections.columns.tolist()[:-1], y=co_infections.columns.tolist()[1:][::-1], colorscale='Viridis'))
    fig.update_layout({"paper_bgcolor": "rgba(0, 0, 0, 0)","plot_bgcolor": "rgba(0, 0, 0, 0)"})
    fig = fig.update_traces(text=co_infections.values[1:,:-1][::-1].astype(str), texttemplate="%{text}", hovertemplate=None, showscale=False)

    co_inf_mat = plot(fig, output_type="div")

    #df=df[~df['strain_name'].isin(['Influenza A', 'Influenza B'])]
    #co_infections=co_infections.drop(index = ['Influenza A', ], columns = ['Influenza A', ])

    if not os.path.isdir('dashboard/static/mixed'):
        os.mkdir('dashboard/static/mixed')
        ref_seqs = pd.read_csv('dashboard/static/reference_seqs.csv')
        df = df.merge(ref_seqs, how='left', left_on = 'substrain_name', right_on='name')
        df['pileup_id'] = df['pseudonymized_id'] + '_' + df['accessionID']

        coverage = pd.read_csv('dashboard/static/pileup_coinfections.csv')
        cols = [col for col in coverage.columns if col in df['pileup_id'].values]
        coverage = coverage[cols]

        for i,j in np.argwhere(np.array(co_infections) > 0):
            if j>i:
                print(co_infections.index[i] + ' + ' +co_infections.index[j])

                samples = mat[(mat[co_infections.index[i]]==1) & (mat[co_infections.index[j]]==1)].index

                if((co_infections.index[i]!='Influenza A') and (co_infections.index[j]!='Influenza A')):
                    fig, (ax1, ax2, ax3, ax4) = plt.subplots(
                        4, 1,  figsize=(20, 20), sharex=False, gridspec_kw={"height_ratios": [4, 1, 4, 1]}
                    )
                    
                    
                    samples_i = df[df['pseudonymized_id'].isin(samples) & (df['strain_name']==co_infections.index[i])]['pileup_id'].tolist()
                    samples_j = df[df['pseudonymized_id'].isin(samples) & (df['strain_name']==co_infections.index[j])]['pileup_id'].tolist()


                    match co_infections.index[i]:
                        case 'Rhino- / Enterovirus':
                            strain_name = 'Enterovirus'
                        case 'RSV - A/B':
                            strain_name = 'RSV'
                        case _:
                            strain_name= co_infections.index[i].replace(' ' ,'_')

                    cols=[s for s in coverage.columns if s in samples_i] 
                    coverage_strain=coverage[cols].fillna(0)
                    coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                    pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                    

                    ax1.plot(pileup['idx'], pileup['mean'])
                    ax1.fill_between(pileup['idx'], pileup['ci_lower'], pileup['ci_upper'], color='b', alpha=.15)
                    ax1.set_ylim(bottom=0)
                    ax1.set_ylabel("Sequencing Depth", fontsize = 10)
                    ax1.set_title(co_infections.index[i],fontsize = 15, loc='left' ,pad=20)
                    ax1.set_yscale("log")
                    ax1.set_ylim(ymin=1)
                    ax1.get_xaxis().set_visible(False)
                    #ax1.axhline(y = 5, color = 'b', linestyle = '--', label = "DP5")
                    ax1.axhline(y = 10, color = 'b', linestyle = '-', label = "DP10")
                    #ax1.axhline(y = 20, color = 'b', linestyle = ':', label = "DP20")
                    ax1.legend(bbox_to_anchor = (1, 1), loc = 'upper left')

                    graphic_record = BiopythonTranslator().translate_record('dashboard/static/annotations/'+strain_name+'.gb')
                    graphic_record.plot(ax=ax2, strand_in_label_threshold=4, with_ruler=True)
                    ax2.set_xlabel("Position")
                    ax2.get_yaxis().set_visible(False)


                    match co_infections.index[j]:
                            case 'Rhino- / Enterovirus':
                                strain_name = 'Enterovirus'
                            case 'RSV - A/B':
                                strain_name = 'RSV'
                            case _:
                                strain_name= co_infections.index[j].replace(' ' ,'_')

                    cols=[s for s in coverage.columns if s in samples_j] 
                    coverage_strain=coverage[cols].fillna(0)
                    coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                    pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                    

                    ax3.plot(pileup['idx'], pileup['mean'])
                    ax3.fill_between(pileup['idx'], pileup['ci_lower'], pileup['ci_upper'], color='b', alpha=.15)
                    ax3.set_ylim(bottom=0)
                    ax3.set_ylabel("Sequencing Depth", fontsize = 10)
                    ax3.set_title(co_infections.index[j],fontsize = 15, loc='left' ,pad=20)
                    ax3.set_yscale("log")
                    ax3.set_ylim(ymin=1)
                    ax3.get_xaxis().set_visible(False)
                    #, color = 'b', linestyle = '--', label = "DP5")
                    ax3.axhline(y = 10, color = 'b', linestyle = '-', label = "DP10")
                    #ax3.axhline(y = 20, color = 'b', linestyle = ':', label = "DP20")
                    ax3.legend(bbox_to_anchor = (1, 1), loc = 'upper left')

                    graphic_record = BiopythonTranslator().translate_record('dashboard/static/annotations/'+strain_name+'.gb')
                    graphic_record.plot(ax=ax4, strand_in_label_threshold=4, with_ruler=True)
                    ax4.set_xlabel("Position")
                    ax4.get_yaxis().set_visible(False)

                else:


                    if(co_infections.index[i]=='Influenza A'):
                        second = co_infections.index[j]
                    else:
                        second = co_infections.index[i]

                    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(
                        6, 1,  figsize=(20, 20), sharex=False, gridspec_kw={"height_ratios": [4, 1, 4, 1, 4, 1]}
                    )
                    
                    samples_HA = [s + '_KU509703.1' for s in df[df['pseudonymized_id'].isin(samples)]['pseudonymized_id']] + [s + '_KJ942616.1' for s in df[df['pseudonymized_id'].isin(samples)]['pseudonymized_id']]
                    samples_NA = [s + '_KU509705.1' for s in df[df['pseudonymized_id'].isin(samples)]['pseudonymized_id']] + [s + '_KJ942618.1' for s in df[df['pseudonymized_id'].isin(samples)]['pseudonymized_id']]
                    samples_second = df[df['pseudonymized_id'].isin(samples) & (df['strain_name']==second)]['pileup_id'].tolist()

                    cols=[s for s in coverage.columns if s in samples_HA] 
                    coverage_strain=coverage[cols].fillna(0)
                    coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                    pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                    

                    ax1.plot(pileup['idx'], pileup['mean'])
                    ax1.fill_between(pileup['idx'], pileup['ci_lower'], pileup['ci_upper'], color='b', alpha=.15)
                    ax1.set_ylim(bottom=0)
                    ax1.set_ylabel("Sequencing Depth", fontsize = 10)
                    ax1.set_title("Influenza A - HA",fontsize = 15, loc='left' ,pad=20)
                    ax1.set_yscale("log")
                    ax1.set_ylim(ymin=1)
                    ax1.get_xaxis().set_visible(False)
                    #, color = 'b', linestyle = '--', label = "DP5")
                    ax1.axhline(y = 10, color = 'b', linestyle = '-', label = "DP10")
                    #ax1.axhline(y = 20, color = 'b', linestyle = ':', label = "DP20")
                    ax1.legend(bbox_to_anchor = (1, 1), loc = 'upper left')

                    graphic_record = BiopythonTranslator().translate_record('dashboard/static/annotations/Influenza_A_HA.gb')
                    graphic_record.plot(ax=ax2, strand_in_label_threshold=4, with_ruler=True)
                    ax2.set_xlabel("Position")
                    ax2.get_yaxis().set_visible(False)

                    cols=[s for s in coverage.columns if s in samples_NA] 
                    coverage_strain=coverage[cols].fillna(0)
                    coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                    pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                    

                    ax3.plot(pileup['idx'], pileup['mean'])
                    ax3.fill_between(pileup['idx'], pileup['ci_lower'], pileup['ci_upper'], color='b', alpha=.15)
                    ax3.set_ylim(bottom=0)
                    ax3.set_ylabel("Sequencing Depth", fontsize = 10)
                    ax3.set_title("Influenza A - NA",fontsize = 15, loc='left' ,pad=20)
                    ax3.set_yscale("log")
                    ax3.set_ylim(ymin=1)
                    ax3.get_xaxis().set_visible(False)
                    #ax3.axhline(y = 5, color = 'b', linestyle = '--', label = "DP5")
                    ax3.axhline(y = 10, color = 'b', linestyle = '-', label = "DP10")
                    #ax3.axhline(y = 20, color = 'b', linestyle = ':', label = "DP20")
                    ax3.legend(bbox_to_anchor = (1, 1), loc = 'upper left')

                    graphic_record = BiopythonTranslator().translate_record('dashboard/static/annotations/Influenza_A_NA.gb')
                    graphic_record.plot(ax=ax4, strand_in_label_threshold=4, with_ruler=True)
                    ax4.set_xlabel("Position")
                    ax4.get_yaxis().set_visible(False)

                    match second:
                            case 'Rhino- / Enterovirus':
                                strain_name = 'Enterovirus'
                            case 'RSV - A/B':
                                strain_name = 'RSV'
                            case _:
                                strain_name= co_infections.index[j].replace(' ' ,'_')

                    cols=[s for s in coverage.columns if s in samples_second] 
                    coverage_strain=coverage[cols].fillna(0)
                    coverage_strain=coverage_strain[:np.where(coverage_strain.sum(axis=1))[0][-1]]
                    pileup= pd.DataFrame({'idx': coverage_strain.index+1, 'mean': coverage_strain.mean(axis=1), 'ci_lower': coverage_strain.mean(axis=1) - 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns)), 'ci_upper': coverage_strain.mean(axis=1) + 1.96 *coverage_strain.std(axis=1)/ np.sqrt(len(coverage_strain.columns))})
                    

                    ax5.plot(pileup['idx'], pileup['mean'])
                    ax5.fill_between(pileup['idx'], pileup['ci_lower'], pileup['ci_upper'], color='b', alpha=.15)
                    ax5.set_ylim(bottom=0)
                    ax5.set_ylabel("Sequencing Depth", fontsize = 10)
                    ax5.set_title(second,fontsize = 15, loc='left' ,pad=20)
                    ax5.set_yscale("log")
                    ax5.set_ylim(ymin=1)
                    ax5.get_xaxis().set_visible(False)
                    #ax5.axhline(y = 5, color = 'b', linestyle = '--', label = "DP5")
                    ax5.axhline(y = 10, color = 'b', linestyle = '-', label = "DP10")
                    #ax5.axhline(y = 20, color = 'b', linestyle = ':', label = "DP20")
                    ax5.legend(bbox_to_anchor = (1, 1), loc = 'upper left')

                    graphic_record = BiopythonTranslator().translate_record('dashboard/static/annotations/'+strain_name+'.gb')
                    graphic_record.plot(ax=ax6, strand_in_label_threshold=4, with_ruler=True)
                    ax6.set_xlabel("Position")
                    ax6.get_yaxis().set_visible(False)

                fig.savefig('dashboard/static/mixed/pileup_'+co_infections.index[i].replace(' ', '_').replace('/', '_')+'_'+co_infections.index[j].replace(' ', '_').replace('/', '_')+'.png', bbox_inches='tight')
             
    #pairs=[{'title': co_infections.index[i] + ' & '+co_infections.index[j], 'plot': 'mixed/pileup_' + str(i) + '_' + str(j) +'.png'} for i,j in np.argwhere(np.array(co_infections) > 0) if j>i]
    #pairs = [co_infections.index[[i,j]].tolist() for i,j in np.argwhere(np.array(co_infections) > 0) if j>i]
    pairs=[{'first': co_infections.index[i], 'second': co_infections.index[j], 'first_clean': co_infections.index[i].replace(' ', '_').replace('/', '_'), 'second_clean': co_infections.index[j].replace(' ', '_').replace('/', '_')} for i,j in np.argwhere(np.array(co_infections) > 0) if j>i]

    return render(request, "mixed.html", {'co_inf_mat': co_inf_mat, 'strains': co_infections.index.tolist(), 'pairs': pairs, 'no_co_infections': no_co_infections }) 