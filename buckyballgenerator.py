


def generator():
    indx = [[5, 6, 2],
            [1, 7, 3],
            [2, 8, 4],
            [3, 9, 5],
            [4, 10, 1],
            [1, 20, 11],
            [2, 12, 13],
            [3, 14, 15],
            [4, 16, 17],
            [5, 18, 19],
            [6, 21, 12],
            [7, 11, 22],
            [7, 23, 14],
            [8, 13, 24],
            [8, 25, 16],
            [9, 15, 26],
            [9, 27, 18],
            [10, 17, 28],
            [10, 29, 20],
            [6, 19, 30],
            [11, 30, 31],
            [12, 32, 23],
            [13, 22, 33],
            [14, 34, 25],
            [15, 24, 35],
            [16, 36, 27],
            [17, 26, 37],
            [18, 38, 29],
            [19, 28, 39],
            [20, 40, 21],
            [21, 41, 32],
            [22, 31, 42],
            [23, 34, 43],
            [24, 33, 44],
            [25, 45, 36],
            [26, 46, 35],
            [27, 47, 38],
            [28, 48, 37],
            [29, 49, 40],
            [30, 50, 39],
            [31, 50, 51],
            [32, 51, 43],
            [33, 42, 52],
            [34, 52, 45],
            [35, 44, 53],
            [36, 53, 47],
            [37, 46, 54],
            [38, 54, 49],
            [39, 48, 55],
            [40, 55, 41],
            [41, 42, 56],
            [43, 44, 57],
            [45, 58, 46],
            [47, 48, 59],
            [49, 50, 60],
            [60, 51, 57],
            [56, 52, 58],
            [57, 53, 59],
            [54, 58, 60],
            [55, 59, 56]]
    for i in range(3):
        for j in range(60):
            indx[j][i] -= 1

    structure = {}
    for i in range(60):
        neigh = tuple(indx[i])
        node = i
        structure[node] = neigh
    edges = [(0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7), (3, 8), (3, 4), (4, 9), (5, 10), (5, 19), (6, 11), (6, 12), (7, 13), (7, 14), (8, 16), (8, 15), (9, 17), (9, 18), (10, 11), (10, 20), (11, 21), (12, 22), (12, 13), (13, 23), (14, 24), (14, 15), (15, 25), (16, 17), (16, 26), (17, 27), (18, 19), (18, 28), (19, 29), (20, 29), (20, 30), (21, 22), (21, 31), (22, 32), (23, 24), (23, 33), (24, 34), (25, 26), (25, 35), (26, 36),
             (27, 28), (27, 37), (28, 38), (29, 39), (30, 40), (30, 31), (31, 41), (32, 33), (32, 42), (33, 43), (34, 35), (34, 44), (35, 45), (36, 37), (36, 46), (37, 47), (38, 48), (38, 39), (39, 49), (40, 49), (40, 50), (41, 50), (41, 42), (42, 51), (43, 51), (43, 44), (44, 52), (45, 52), (45, 46), (46, 53), (47, 48), (47, 53), (48, 54), (49, 54), (50, 55), (51, 56), (52, 57), (53, 58), (54, 59), (55, 56), (55, 59), (56, 57), (57, 58), (58, 59)]

    return structure, edges