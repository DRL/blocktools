def block_algorithm_b_old(params):
    regionBatchObj, parameterObj = params
    block_idx = 0
    blockObjs = deque()
    block_id = "%s.i%s.b%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx)
    blockObj = BlockObj(block_id, parameterObj.block_length)
    #print("\n[NEW]", blockObj)
    while 1:
        try:
            bedObj = regionBatchObj.bedObjs.popleft()
            #print("[->]", bedObj)
            if blockObj.pair_idxs:
                overlap_pair_idxs = compatible_pair_idxs(blockObj.pair_idxs, bedObj.pair_idxs)
                #print("[OVL]", overlap_pair_idxs)
                if not overlap_pair_idxs:
                    #print("[->]", bedObj)
                    bedObj = regionBatchObj.bedObjs.popleft()
                    #print("\n[NEW]", blockObj)
                    blockObj = BlockObj(block_id, parameterObj.block_length)
                else:
                    bedObj.pair_idxs = overlap_pair_idxs
            #if len(bedObj.pair_idxs) >= 2:
            bedObj = blockObj.add_bedObj(bedObj, parameterObj)
            if blockObj.blocked:
                #print("[B!]", blockObj)
                blockObjs.append(blockObj)
                block_idx += 1
                block_id = "%s.i%s.b%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx)
                blockObj = BlockObj(block_id, parameterObj.block_length)
            if bedObj:
                #print("[<-]", bedObj)
                regionBatchObj.bedObjs.appendleft(bedObj)
        except IndexError:
            break
    return blockObjs

def compatible_pair_idxs(block_pair_idxs, bed_pair_idxs):
    '''block_score'''

    difference = bed_pair_idxs.difference(block_pair_idxs)
    intersection = block_pair_idxs.intersection(bed_pair_idxs)
    if len(intersection) >= len(difference):
        return intersection
    else:
        return {}

def block_algorithm_b_2(params):

    '''
    d_i = (bedObj.pair_idxs / parameterObj.pairs_count) * (bedObj.length / parameterObj.block_length)
    '''
    regionBatchObj, parameterObj = params
    block_idx = 0
    blockObjs = deque()
    _blockObjs = []
    #print("\n####", regionBatchObj.contig_id, regionBatchObj.idx)
    #for bedObj in regionBatchObj.bedObjs:
    #    print(bedObj)ew
    while 1:
        try:
            bedObj = regionBatchObj.bedObjs.popleft()
            #print(bedObj)
            if not _blockObjs:
                block_idx += 1
                block_id = "%s.i%s.b%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx)
            else:
                _blockObjs = [_blockObj for _blockObj in _blockObjs if _blockObj.score > 0.0]
            blockObj = BlockObj(block_id, parameterObj.block_length)
            _blockObjs.append(blockObj)
            for idx, _blockObj in enumerate(sorted(_blockObjs, key=lambda i: (i.score), reverse=True)):
                remainder_bedObj = _blockObj.add_bedObj(bedObj, parameterObj)
                #print("=> \t %s" % _blockObj)
                if blockObj.score > 0.0:
                    if blockObj.needed == 0:
                        blockObjs.append(_blockObj)
                        _blockObjs = []
                        #print("# final: ", "\n%s" % _blockObj)
                        if remainder_bedObj:
                            regionBatchObj.bedObjs.appendleft(remainder_bedObj)
        except IndexError:
            break
    return blockObjs
    
def block_algorithm_c(params):
    '''
    - able to jump gaps in intervals
    - able to jump intervals

    a) bedObj_score : score of bedObj if a new block is started
    b) list of single bedObjs that have been pop'ped, to populate blockObj candidates
    c) list of blockObj_candidates
        - composed of blockObjs with bedObj
        - blockObj
    condition for survival in c:
        - c = [blockObj for blockObj in c if blockObj.score > 0]
        - c = [blockObj for blockObj in c if blockObj.score >= bedObj_score]

    bedObj = bedObjs.popleft()
    if c:
        c_copy = c.deep_copy()
        for _c in c:
            _c.add_bedObj(bedObj)
        c.append(blockObj)
        c.extend(c.copy)
    else:
        blockObj.append(bedObj)
        c.append(blockObj)
    c = [blockObj for blockObj in c if blockObj.score > 0]    
    
    candidate for adding bedObjs

    i = 1: 
        seeds = [ ]  
        candidates = [1] if blockObj.score > bedObj.score
        pass_on = [1]
        finished = []
        remainders = []
    i = 2: 
        seeds = [1]
        candidates = [12, 2] if blockObj.score > bedObj.score
        pass_on = [1, 12, 2]
        finished = []
        remainders = []
    i = 3: 
        seeds = [1, 12, 2] ;
        candidates = [13, 123, 23, 3] if blockObj.score > bedObj.score
        pass_on = [1, 12, 2, 13, 123, 23, 3]
        finished = []
        remainders = []
    i = 4: 
        seeds = [1, 12, 2, 13, 123, 23, 3] 
        candidates = [14, 124, 24, 134, 1234, 234, 4] if blockObj.score > bedObj.score 
        finished = [1234, 234]
        remainders = [4", 4""]
    i = 5: 
        seeds = [1, 12, 2, 13, 123, 23, 3, 14, 124, 24, 134, 1234, 234, 4] 
        candidates = [15, 125, 25, 1235, 235, 35, 145, 1245, 245, 1345, 12345, 2345, 45, 5] if blockObj.score > bedObj.score
        finished = [1234, 234, 235, 2345, 45]
        remainders = [4", 4"", 5", 5"", 5"""]

    Ways of removing seeds:
        - span violation
        - span_rate violation 
        - score > bedObj.score  

    way of limiting finished:
        - max_finished?
        - min_score?

    
    '''

    regionBatchObj, parameterObj = params
    block_idx = 0
    blockObjs = deque()
    block_id = "%s.i%s.b%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx)
    blockObj = BlockObj(block_id, parameterObj.block_length)
    candidates = []
    finished = []
    remainders = []
    #print("\n", regionBatchObj.contig_id, regionBatchObj.idx)
    #print("\n[NEW]", blockObj)
    while 1:
        try:
            bedObj = regionBatchObj.bedObjs.popleft()
            # blockObj.needed of current blockObj determines bedObj_score
            #bedObj_score = (len(bedObj.pair_idxs) / parameterObj.pairs_count) * (min(bedObj.length, blockObj.needed) / parameterObj.block_length)
            # bedObj is added to blockObj, remainder bedObj is returned
            #print("[BED]", bedObj)
            new_candidates = []
            if candidates:
                #print("candidates ...")
                while candidates:
                    blockObj = candidates.pop()
                    old_blockObj = copy.deepcopy(blockObj)
                    new_candidates.append(old_blockObj)
                    remainder_bedObj = blockObj.add_bedObj(bedObj, parameterObj)
                    if blockObj.needed == 0 and blockObj.score > 0.0:
                        #print("\t[F]", blockObj)
                        finished.append(blockObj)
                        if remainder_bedObj:
                            remainders.append(remainder_bedObj)
                    else:
                        #print("\t[C]", blockObj)
                        new_candidates.append(blockObj)
            else:
                #print("starting new block ...")
                blockObj = BlockObj(block_id, parameterObj.block_length)
                #print("\t[N]", blockObj)
                remainder_bedObj = blockObj.add_bedObj(bedObj, parameterObj)
                if blockObj.needed == 0 and blockObj.score:
                    #print("\t[F]", blockObj)
                    finished.append(blockObj)
                    if remainder_bedObj:
                        remainders.append(remainder_bedObj)
                else:
                    new_candidates.append(blockObj)
                    #print("\t[C]", blockObj)
            if finished:
                #print("finished ...")
                max_value, max_idx = 0.0, 0 
                for idx, blockObj in enumerate(finished):
                    #print("\t[F]: ", blockObj)
                    if blockObj.score >= max_value:
                        max_value = blockObj.score
                        max_idx = idx
                    #print("[FB]", idx, finished[idx])
                #print("\t[*]: ", blockObj)
                blockObjs.append(finished[max_idx])
                #print("[*]", idx, finished[idx])
                finished = []
                candidates = []
                block_idx += 1
                blockObj = BlockObj(block_id, parameterObj.block_length)
                try:
                    if remainders[max_idx]:
                        #print("\t<-", remainder_bedObj)
                        regionBatchObj.bedObjs.appendleft(remainders[max_idx])
                    remainders = []
                except IndexError:
                    pass
            #candidates = [blockObj for blockObj in new_candidates if blockObj.score > 0.0]
        except IndexError:
            break
        if finished:
                #print("finished ...")
                max_value, max_idx = 0.0, 0 
                for idx, blockObj in enumerate(finished):
                    #print("\t[F]: ", blockObj)
                    if blockObj.score >= max_value:
                        max_value = blockObj.score
                        max_idx = idx
                    #print("[FB]", idx, finished[idx])
                #print("\t[*]: ", blockObj)
                blockObjs.append(finished[max_idx])
                #print("[*]", idx, finished[idx])
                #try:
                #    if remainders[max_idx]:
                        #print("\t<-", remainder_bedObj)
                #        regionBatchObj.bedObjs.appendleft(remainders[max_idx])
                #except IndexError:
                #    pass
                #finished = []
                #remainders = []
                #candidates = []
                #block_idx += 1
                #blockObj = BlockObj(block_id, parameterObj.block_length)
    return blockObjs

def findMaxAverage(arr, n, k): 
    # Check if 'k' is valid 
    if (k > n): 
        return -1
  
    # Compute sum of first 'k' elements 
    _sum = sum([value for value in range(0, k)])  
  
    max_sum = _sum
    max_end = k - 1
  
    # Compute _sum of remaining subarrays 
    for i in range(k, len(arr)): 
        _sum = _sum + arr[i] - arr[i - k]  
        if (_sum > max_sum): 
            max_sum = _sum
            max_end = i  
          
    # Return starting index 
    return max_end - k + 1

#def generate_blockObj(parameterObj, bedObjs):


            #for bedObj in bedObjs:
            #length += bedObj.length
            #if length > parameterObj.block_length:
        #span = bedObj
        #if sum([bedObj.length for bedObj in bedObjs])     
        #span = [bedObj.length bedObjs
        #if not blockObj.pair_idxs: # fresh block, no debate
            #bedObj = blockObj.add_bedObj(bedObj, parameterObj)
        #else:
            #required_length = parameterObj.block_length - blockObj.length
            #B_i = min(bedObj.length, required_length)
            #d_1 = (len(bedObj.pair_idxs) / parameterObj.pairs_count) * (B_i / parameterObj.block_length)
            #d_2 = 
        #if not blockObj.pair_idxs: # fresh block, no debate
            #bedObj = blockObj.add_bedObj(bedObj, parameterObj)
            #print(blockObj)
        #else:
            #required_length = parameterObj.block_length - blockObj.length
            #P_n = blockObj.pair_idxs.intersection(bedObj.pair_idxs)
            #print("[needed]", required_length, "intersection", P_n, len(P_n))
            #B_n = blockObj.length
            #d_n = (len(P_n) / parameterObj.pairs_count) * (B_n / parameterObj.block_length)
            #B_i = min(bedObj.length, required_length)
            #d_i = (len(bedObj.pair_idxs) / parameterObj.pairs_count) * (B_i / parameterObj.block_length)
            #if d_n > d_i:
                #print("[d_n]", d_n, ">", "[d_i]", d_i)
                #bedObj = blockObj.add_bedObj(bedObj, parameterObj)
                #print(blockObj)
            #else:
                #print("[d_n]", d_n, "<", "[d_i]", d_i)
                #bedObj = regionBatchObj.bedObjs.popleft()
                #print("[->]", bedObj)
            #print("\n[NEW]", blockObj)
            #blockObj = BlockObj(block_id, parameterObj.block_length)            
        #if blockObj.blocked:  # blocked when full
            #print("[@]", blockObj)
            #blockObjs.append(blockObj)
            #block_idx += 1
            #block_id = "%s.i%s.b%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx)
            #blockObj = BlockObj(block_id, parameterObj.block_length)
        #if bedObj:
            #print("[<-]", bedObj)
    #regionBatchObj.bedObjs.appendleft(bedObj)