import os
import pandas as pd
import numpy as np
from collections import defaultdict
hifi_file = '/Users/xrk.study/Desktop/bgi_graduate/MF2.verkko.mat.map2t2t.paf'

info_to_extract = ['q_id'  , 'q_length' , 'q_st' , 'q_stp' ,'char' ,  'r_name' , \
                           'r_length' , 'r_st' , 'r_stp' ,'num_match' , 'num_bases' , 'quality' , 'tp' , 'mm' ,
                   'gn' , 'go' , 'cg']
output_agp = '/Users/xrk.study/Desktop/MF2.verkko.mat.map2t2t.agp'
os.system('touch {}'.format(output_agp) )

def find_insertion_sets(set_list ) : # function to merge sets
    merge_sets = []
    if len(set_list) == 0 :
        print('input sets empty')
        return
    else :
        a1 , b1 = set_list[0][0] , set_list[0][1]
        for i in range(1 , len(set_list)) :
            a2 , b2 = set_list[i][0] , set_list[i][1]
            # 只要 a2 < b1 就有重叠
            if b1 >= a2 :
                b1 = max(b1 , b2)
            else :
                merge_sets.append([a1 , b1])
                a1 , b1 = a2 , b2
    merge_sets.append([a1 , b1])
    return merge_sets

def get_subcont_char(list , in_dict) :

    return


def clean_short_sets(dict , merge_sets_min ) :   # filter to clean short alignment
    clean_dict = {}
    for k , v in dict.items() :
        merge_sets_length = v['q_sets'][0][1] - v['q_sets'][0][0]
        if merge_sets_length <= merge_sets_min:
            continue
        else:
            clean_dict[k] = v

    return clean_dict



def extra_dict_to_list(in_dict ,  key ) :
    ou_list = []
    for k, v in in_dict.items() :
        ou_list.append(v[key])

    return ou_list


contig_length_min = 100 # contigs >= 10 kb
# initialize
query_id = ''
ref_chr = ''
all_query_ids = []
all_chrs = []
# all alignment
align_sum = {}
ouptut_sum = {}


def sort_by_ref(in_dict):
    out_dict = {}
    ref_coords = sorted([v['r_merge_sets'][0] for k, v in in_dict.items()] )
    sort_ref_coords = []
    for rc in ref_coords:
        sort_ref_coords.append([v for k, v in in_dict.items() if v['r_merge_sets'][0] == rc])

    for i in range(len(sort_ref_coords)) :
        out_dict['set_{}'.format(str(i))] = sort_ref_coords[i]

    return out_dict


def append_break_list(in_list , index) :

    ou_list = []
    ou_list.append(in_list[index])
    ou_list.append(in_list[index + 1 ])

    return ou_list



with open(hifi_file , 'r+') as f :
    for line in f :
        lineinfo = line.strip().split()
        lineid = lineinfo[0]
        if lineid not in align_sum.keys() :
            id_dict = {}
            for i in range(len(info_to_extract)) :
                id_dict[info_to_extract[i]] = [lineinfo[i] ]
        elif lineid in align_sum.keys() :
            for i in range(len(info_to_extract)) :
                id_dict[info_to_extract[i]].append(lineinfo[i])

        align_sum[lineid] = id_dict

all_agp_list = []

# ###### params could be given use command line

# print(align_sum)
#

min_merge_align = 10000  # do not consider alignment < 10 kb
find_chimeric = True  # find chimeric contigs
close_cutoff = 0  # define close cutoff
min_coverage_rate = 0.8 # min coverage rate
chimeric_coverage_range = 0.05
chimeric_set_distance = 100000  # 100 kb means chimeric contigs

for k, v in align_sum.items() :
    # 判断 有结果 align
    v['q_st']  =[int(x) for x in v['q_st']]
    v['q_stp'] =[int(x) for x in v['q_stp']]
    v['r_st'] =  [int(x) for x in v['r_st']]
    v['r_stp'] = [int(x) for x in v['r_stp']]
    v['q_coords'] = [[v['q_st'][i] , v['q_stp'][i]] for i in range(len(v['q_st']))]
    v['r_coords'] = [[v['r_st'][i], v['r_stp'][i]]  for i in range(len(v['q_st']))]

    if len(v['q_id']) == 1 :
        contig_name = v['q_id'][0]
        q_length = v['q_length'][0]
        char = v['char'][0] # char
        chr =  v['r_name'][0]
        coverage = (int(v['q_stp'][0] ) - int(v['q_st'][0])) / int(q_length)
        if coverage > min_coverage_rate :


            agp_info = contig_name + '\t1\t{}\t1\tW\t{}_subcontig1\t1\t{}\t{}\t{}'.format(q_length,
                                                                       contig_name,
                                                                       q_length,
                                                                       char , chr)
            all_agp_list.append(agp_info)
        else :
            continue


    if len(v['q_id']) > 1 :
        align_refs = list(set(v['r_name']) )
        if len(align_refs) > 1 :   #
            chr_info_dict = {}
            contig_name = v['q_id'][0]
            q_length = int(v['q_length'][0] )

            all_chrs_info = []

            chimeric_ref_chr = []
            chimeric_qsets = []
            chimeric_contig_length = []

            for chr in align_refs :
                chr_dict = {}
                q_coords = [[int(v['q_st'][i]), int(v['q_stp'][i])] for i in range(len(v['q_st'])) if v['r_name'][i] == chr] # chr
                r_coords = [[int(v['r_st'][i]), int(v['r_stp'][i])] for i in range(len(v['q_st'])) if v['r_name'][i] == chr]
                char = [v['char'][i] for i in range(len(v['q_st'])) if v['r_name'][i] == chr ]   ## char
                identity = [int(v['num_match'][i]) / int(v['num_bases'][i]) for i in range(len(v['q_st'])) if v['r_name'][i] == chr ]  # identity
                total_length = 0
                merge_qsets = find_insertion_sets(sorted(q_coords))  # merge sets
                for ms in merge_qsets :
                    total_length += ms[1] - ms[0]
                total_coverage_rate = total_length / q_length

                chr_dict['contig_name'] = contig_name
                chr_dict['q_coords'] = q_coords
                chr_dict['q_merge_coords'] = find_insertion_sets(sorted(q_coords ) )
                chr_dict['q_length'] = q_length
                chr_dict['r_coords'] = r_coords
                chr_dict['r_merge_coords'] = find_insertion_sets(sorted(r_coords) )
                chr_dict['identity'] = identity
                chr_dict['coverage'] = total_coverage_rate
                chr_dict['ref_chr'] = chr
                chr_dict['char'] = char


                all_chrs_info.append(chr_dict)
                if total_coverage_rate > min_coverage_rate :
                    chimeric_ref_chr.append(chr)
                    chimeric_qsets.append(chr_dict['q_merge_coords'])
                    chimeric_contig_length.append((len(chr_dict['q_merge_coords']) , len(chr_dict['r_merge_coords'])))

                else :
                    continue


            if len(chimeric_ref_chr) == 0 :
                continue

            if len(chimeric_ref_chr ) == 1 :
                agp_info = contig_name + '\t1\t{}\t1\tW\t{}_subcontig1\t1\t{}\t{}\t{}'.format(q_length,
                                                                                           contig_name,
                                                                                           q_length,
                                                                                          max(char,key=char.count) ,chr)
                all_agp_list.append(agp_info)    # append agp list


            if len(chimeric_ref_chr) == 2  :
                if len(chimeric_contig_length) == 2 :
                    if len(chimeric_contig_length) == 2 and list(set(chimeric_contig_length)) == [(1, 1)] :
                        print(all_chrs_info)
                        q_merge_coords = [x['q_merge_coords'][0] for x in all_chrs_info ]   # chrs info
                        sorted_q_coords = sorted(q_merge_coords)  # sorted
                        char0 = [x['char'] for x in all_chrs_info if x['q_merge_coords'] == sorted_q_coords[0]]
                        char1 = [x['char'] for x in all_chrs_info if x['q_merge_coords'] == sorted_q_coords[1]]
                        chr0 = [x['ref_chr'] for x in all_chrs_info if x['q_merge_coords'] == sorted_q_coords[0]]
                        chr1 = [x['ref_chr'] for x in all_chrs_info if x['q_merge_coords'] == sorted_q_coords[1]]

                        agp_info1 = contig_name + '\t1\t{}\t1\tW\t{}_subcontig1\t1\t{}\t{}\t{}'.format(sorted_q_coords[0][1] ,
                                                                                                   contig_name  , sorted_q_coords[0][1] ,
                                                                                                   max(char0 , key=char0.count) , chr0)
                        agp_info2 = contig_name + '\t{}\t{}\t2\tW\t{}_subcontig2\t1\t{}\t{}\t{}'.format(sorted_q_coords[0][1] ,
                                                                                                   q_length , contig_name , (int(q_length) - int(sorted_q_coords[0][1] )  )
                                                                                                   , max(char1 , key=char1.count) , chr1

                                                                                                   )
                        all_agp_list.append(agp_info1)
                        all_agp_list.append(agp_info2)

                      # 这种是需要打断的
                    else :
                        continue

                else :
                    continue
            else :
                continue

               # 只讨论 两种情况下的 打断



        ## shen zhen # #
        if len(align_refs ) == 1 :  # 仅对应知align 到 一个 chr 的 情况
            q_length = int(v['q_length'][0])   # contig length

            if q_length <= contig_length_min :
                continue
            else :
                q_coords =[ [int(v['q_st'][i]) , int(v['q_stp'][i])] for i in range(len(v['q_st'])) ]
                q_coords_merge = find_insertion_sets(sorted(q_coords))  # sorted by query coords
            # 通过 dict 找到所有的 信息
                q_all_sets_info = {}
                for q_m_set in q_coords_merge :

                    merge_info = {}
                    start_coord, stop_coord = q_m_set[0], q_m_set[1]
                    q_sets_to_merge = [x for x in q_coords if x[0] >= start_coord and x[1] <= stop_coord  ]

                    merge_info['q_sets'] = q_sets_to_merge

                    '''
                    start_index = [sorted(q_coords).index(x) for x in sorted(q_coords) if x[0] == start_coord]
                    stop_index = [sorted(q_coords).index(x) for x in sorted(q_coords) if x[1] == stop_coord]
                    for i in range(start_index[0] , stop_index[0] + 1 ) :
                        q_sets_to_merge.append(sorted(q_coords)[i])
                    merge_info['q_sets'] = q_sets_to_merge 
                    '''
               #  merge_info['q_sets_merge'] = find_insertion_sets(q_sets_to_merge)

                    r_sets_to_merge = []
                    coverage = []
                    identity = []
                    char = []
                    contig_name = []
                    chr = ''
                    for s in q_sets_to_merge :
                        index = v['q_coords'].index(s)
                        r_sets_to_merge.append(v['r_coords'][index])
                        coverage.append(int(v['num_bases'][index]) / int( q_length) )
                        identity.append(int(v['num_match'][index]) / int(v['num_bases'][index]) )
                        char.append(v['char'][index])
                        contig_name.append(v['q_id'][index])
                        chr = v['r_name'][0]

                    merge_info['contig_name'] = contig_name[0]

                    merge_info['coverage'], merge_info['identity'] = coverage , identity
                    # char 为 这些 set 中 出现 次数 最多 的
                    merge_info['char'] = max(char, key=char.count )

                    # q_contig name

                    merge_info['r_sets'] = r_sets_to_merge
                    # merge ref set
                    merge_info['r_sets_merge'] = find_insertion_sets(sorted(r_sets_to_merge))
                    merge_info['q_length'] = q_length
                    merge_info['ref_chr'] = chr

                    q_all_sets_info['set_{}'.format(q_coords_merge.index(q_m_set))] = merge_info


                   #  print(q_all_sets_info)

                '''
                    if merge_info['r_sets'][0][0] > merge_info['r_sets'][0][1] :
                        continue 
                    else : 
                        print(merge_info)
                '''
                q_merge_dict = defaultdict(dict)

                # print(q_all_sets_info)

                for k1 ,  v1 in q_all_sets_info.items() :
                    #  当 这个 sets 仅 对应 ref 上面 一段 set 的 时候
                    if len(v1['r_sets_merge'] )  == 1 :
                        q_merge_sets = find_insertion_sets(sorted(v1['q_sets']))
                        if q_merge_sets[0][1] - q_merge_sets[0][0] <  min_merge_align :
                            continue

                        else :

                            q_merge_dict[k1]['q_merge_sets'] = q_merge_sets
                            q_merge_dict[k1]['coverage'] =(  q_merge_sets[0][1] - q_merge_sets[0][0] )  / q_length
                            q_merge_dict[k1]['identity'] = sum(v1['identity']) / len(v1['identity'] )
                            q_merge_dict[k1]['r_merge_sets'] = v1['r_sets_merge']
                            q_merge_dict[k1]['q_length'] = q_length
                            q_merge_dict[k1]['char'] = v1['char']
                            q_merge_dict[k1]['contig_name'] = v1['contig_name']
                            q_merge_dict[k1]['ref_chr'] = v1['ref_chr']

                    # 当 这个 query 区间 对应两个 及以上的 区间的时候，保留 coverage 最大的 那个
                    else :
                        q_merge_sets = find_insertion_sets(sorted(v1['q_sets']))  # q merge dict
                        r_sets , r_sets_merge = v1['r_sets'] , v1['r_sets_merge']   # r sets merge

                        # 选出 merge align length > min_merge_align 的
                        r_sets_filter_sets = [x for x in r_sets_merge if x[1] - x[0] > min_merge_align]  # 排除 长度小于 10kb 的 alignment
                        if len(r_sets_filter_sets) < 1 :
                            continue
                        else :

                            merge_identity , merge_coverage  = [] , []
                            identity , coverage = v1['identity'] , v1['coverage']
                            for rms in r_sets_filter_sets  :

                                merge_coverage.append((rms[1] - rms[0]) / q_length)
                                merge_rsets = [x for x in r_sets if x[0] >= rms[0] and x[1] <= rms[1]]  # rms 1
                                merge_rsets_identity_list = [identity[r_sets.index(x)] for x in merge_rsets]  # identity
                                merge_rsets_identity = sum(merge_rsets_identity_list) / len(merge_rsets_identity_list)
                                merge_identity.append(merge_rsets_identity)  # append merge rsets

                            max_coverage = max(merge_coverage)     # max coverage rate we consider
                            max_coverage_index = merge_coverage.index(max_coverage)
                            max_coverage_identity = merge_identity[max_coverage_index]
                            max_coverage_rmsets  = r_sets_filter_sets[max_coverage_index]


                            q_merge_dict[k1]['q_merge_sets'] = q_merge_sets
                            q_merge_dict[k1]['coverage'] = max_coverage  ## merge coverage
                            q_merge_dict[k1]['identity'] = max_coverage_identity  ## merge identity
                            q_merge_dict[k1]['r_merge_sets'] =[ max_coverage_rmsets ]
                            q_merge_dict[k1]['q_length'] = q_length
                            q_merge_dict[k1]['char'] = v1['char']
                            q_merge_dict[k1]['contig_name'] = v1['contig_name']
                            q_merge_dict[k1]['ref_chr'] = v1['ref_chr']


                if len(q_merge_dict.keys()) == 0:
                    # 这种事 contig 太短，不考虑
                    continue
                   #  print('{} contigs are too short !'.format(k))
                else :

                    total_cover_length = 0
                    for k , v in q_merge_dict.items() :
                        total_cover_length += v['q_merge_sets'][0][1] - v['q_merge_sets'][0][0]
                    total_coverage_rate = total_cover_length / q_length

                    # min coverage rate 之下的 认为没有 比对上
                    # 这种情况是直接没 比对上，不考虑
                    if total_coverage_rate < min_coverage_rate :
                        continue

                    else :
                        covs = extra_dict_to_list(q_merge_dict, 'coverage')
                        # dair u
                        max_cov = max(covs)
                        if max_cov < 0.1 :
                            continue
                        else :
                            covs_candidate = [x for x in covs if x >= max_cov - chimeric_coverage_range and x <= max_cov]

                            max_candidate_index = covs_candidate.index(max_cov)  # max coverage index
                        #  max2_candidate_index = covs_candidate.index(sorted(covs_candidate)[-2]) # 只是找到 cov 第二的一段
                            if len(covs_candidate) == 1 :   # 如果 没有 相近的 set 或者只有 一个 set 这种时候 不需要 打断
                                # 不需要打断，agp file 这个 contig 只有一段
                                char_list = [xx['char'] for yy, xx in q_merge_dict.items()]
                                agp_char = max(char_list, key=char_list.count)
                                agp_chr = [xx['ref_chr'] for yy , xx in q_merge_dict.items()][0]
                                agp_info = contig_name[0] + '\t1\t{}\t1\tW\t{}_subcontig1\t1\t{}\t{}\t{}'.format(q_length,
                                                                                                       contig_name[0],
                                                                                                       q_length,
                                                                                              agp_char , agp_chr)
                                all_agp_list.append(agp_info)



                            else :

                                sort_ref_merge_dict = sort_by_ref(q_merge_dict)

                           #  print(sort_ref_merge_dict)
                                all_qsets = [v['q_merge_sets'][0] for k ,v in q_merge_dict.items() ]
                                candidate_setid = [k for k , v in sort_ref_merge_dict.items() if v[0]['coverage'] in covs_candidate ]
                                candidate_qsets = [ sort_ref_merge_dict[x][0]['q_merge_sets'][0] for x in candidate_setid ]   # candidate q_merge_sets
                                candidate_rsets = [ sort_ref_merge_dict[x][0]['r_merge_sets'][0] for x in candidate_setid ]

                           #  print(candidate_qsets)
                            # print(candidate_rsets)

                            # 选出了比较 相近的 coverage 对应的 ref query 区间，如果 ref 上面区间差距大就 打断
                           #  print('=====' )
                           #  print('==')
                            # 当 正常的 query 上面这一段
                                distance_ref = []

                                for i in range(len(candidate_rsets) - 1) :
                                    set_d = candidate_rsets[i+1][0] - candidate_rsets[i][1]
                                    distance_ref.append(set_d)

                            # 最大距离的 gap index
                            # print(distance_ref)
                                max_distance_index = distance_ref.index(max(distance_ref) )

                            # chimeric_set_distance = 100000  # 100 kb means chimeric contigs
                            # 当 contig 之间的 distance 小于我们设定的阈值的 时候
                                if max(distance_ref) < chimeric_set_distance :
                                    char_list = [xx['char'] for yy, xx in q_merge_dict.items()]
                                    agp_char = max(char_list, key=char_list.count)
                                    agp_chr = [xx['ref_chr'] for yy, xx in q_merge_dict.items()][0]
                                    agp_info = contig_name[0] + '\t1\t{}\t1\tW\t{}_subcontig1\t1\t{}\t{}\t{}'.format(q_length,
                                                                                                             contig_name[0] ,
                                                                                                             q_length  , agp_char , agp_chr )
                                    all_agp_list.append(agp_info)
                                else :

                               #  porint(sort_ref_merge_dict)
                                    ref_chr = [vv['ref_chr'] for kk , vv in q_merge_dict.items()][0]
                                    candidate_to_break_ref = append_break_list(candidate_rsets , max_distance_index )
                                    candidate_to_break_qry = sorted([candidate_qsets[candidate_rsets.index(x) ] for x in candidate_to_break_ref ] )
                                    subcontig2 = (candidate_to_break_qry[1][0] , candidate_to_break_qry[1][1] )  # subcontig set
                                    subcontig1 = [0 , candidate_to_break_qry[1][0]]   # subcontig2
                                    subcontig1_set = [x for x in all_qsets if x[1] <= subcontig1[1]]
                                    subcontig2_set = [x for x in all_qsets if x[0] >= subcontig2[0]]
                                    subcontig1_set_char = [vv['char'] for kk ,vv in q_merge_dict.items() if vv['q_merge_sets'][0] in subcontig1_set]    # subcontig1_char
                                    subcontig2_set_char = [vv2['char'] for kk2 , vv2 in q_merge_dict.items() if vv2['q_merge_sets'][0] in subcontig2_set]
                                    subcontig1_char = max(subcontig1_set_char , key=subcontig1_set_char.count)
                                    subcontig2_char = max(subcontig2_set_char , key=subcontig2_set_char.count)

                                ## write in cut chimeric contig info
                                    agp_info1 = contig_name[0] + '\t1\t{}\t1\tW\t{}_subcontig1\t1\t{}\t{}\t{}'.format(candidate_to_break_qry[1][0], contig_name[0] ,
                                                                                                          candidate_to_break_qry[1][0] , subcontig1_char ,ref_chr)
                                    agp_info2 = contig_name[0] + '\t{}\t{}\t2\tW\t{}_subcontig2\t1\t{}\t{}\t{}'.format(candidate_to_break_qry[1][0] + 1
                                                                                                               ,q_length ,
                                                                                                                contig_name[
                                                                                                                        0] ,
                                                                                                                q_length  - (candidate_to_break_qry[1][0] + 1 ) ,
                                                                                                               subcontig2_char , ref_chr)
                                    all_agp_list.append(agp_info1)
                                    all_agp_list.append(agp_info2)




with open(output_agp , 'r+') as f:
    for x in all_agp_list :
        f.write(x + '\n')

