"""
Project Name: CSCI 570 Final Project
Author: Luke Nelson, Peter Argo
"""


def string_generator(string: str, jk: list):
    """
    Gererate the stings from the input files
    :param string: input string
    :param jk: list of multipliers 
    :return:
    """
    if len(jk) == 1:
        return string[0:jk[0]+1] + string + string[jk[0]+1::]

    new_string = string[0:jk[0]+1] + string + string[jk[0]+1::]
    return string_generator(new_string, jk[1::])


if __name__ == '__main__':
    s1 = "ACTG"     # "ACACTGACTACTGACTGGTGACTACTGACTGG"
    j = [3,6,1]
    s2 = "TACG"     # "TATTATACGCTATTATACGCGACGCGGACGCG"
    k = [1, 2, 9]
    new_s1 = string_generator(string=s1, jk=j)

    print(new_s1)

