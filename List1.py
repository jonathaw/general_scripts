class List1():
    def __init__(self, items):
        self.list1 = [0] + items

    @property
    def __class__(self):
        return List1

    def __getitem__(self, item):
        """
        :param item: list1 index
        :return: item at list index item(+1)
        >>> l = List1([1, 2, 3])
        >>> l[1]
        1
        >>> l[3]
        3
        """
        return self.list1[item]

    def __setitem__(self, key, value):
        """
        :param key: a postion where to insert value
        :param value: value
        >>> l = List1([1, 2, 3])
        >>> l[2] = 5
        >>> l == List1([1, 5, 3])
        True
        """
        self.list1[key] = value

    def __eq__(self, other):
        """
        :param other: a List1 instance to compare
        :return:True/False
        >>> l1 = List1([1, 2, 3])
        >>> l2 = List1([1, 2, 3])
        >>> l1 == l2
        True
        >>> l3 = List1([1, 4, 3])
        >>> l1 == l3
        False
        """
        for i in range(len(self.list1)):
            if self.list1[i] != other.list1[i]:
                return False
        return True

    def __iter__(self):
        return self.list1[1:].__iter__()

    def __len__(self):
        """
        :return: int length of the list
        >>> l = List1([1, 2, 3])
        >>> len(l)
        3
        """
        return len(self.list1)-1

    def __repr__(self):
        return '<%s>' % ', '.join([str(a) for a in self.list1[1:]])

    def __add__(self, other):
        self.list1 += other

    def __contains__(self, item):
        """
        :param item: an item
        :return: is the item in list
        >>> l = List1([1, 2, 3])
        >>> 1 in l
        True
        >>> 5 in l
        False
        """
        return item in self.list1[1:]

    def append(self, item):
        self.list1.append(item)

    def pop(self, item):
        """
        :param item: item to pos
        :return: the item in item, remove item from list
        >>> l = List1([1, 2, 3])
        >>> l.pop(1)
        1
        >>> l == List1([2, 3])
        True
        """
        return self.list1.pop(item)

    def count(self, item):
        return self.list1[1:].count(item)

    def remove(self, item):
        self.list1.remove(item)

    def index(self, item):
        """
        :param item: item
        :return: index if item
        >>> l = List1([1, 2, 3])
        >>> l.index(2) == 2
        True
        """
        return self.list1.index(item)

    def enumerate(self):
        """
        :return: index, item of List1
        """
        for i in range(1, len(self.list1)):
            yield i, self.list1[i]

    def reverse(self):
        """
        >>> l = List1([1, 2, 3])
        >>> l.reverse()
        >>> l == List1([3, 2, 1])
        True
        """
        self.list1.pop(0)
        self.list1.append(0)
        self.list1.reverse()

    def insert(self, key, value):
        """
        :param key: index
        :param value: value
        :return: List1 with value inserted at key
        >>> l = List1([1, 2, 3])
        >>> l.insert(2, 10)
        >>> l == List1([1, 10, 2, 3])
        True
        """
        self.list1.insert(key, value)

    def sort(self, *kwargs):
        """
        >>> l = List1([1, 2, 4, 3])
        >>> l.sort()
        >>> l == List1([1, 2, 3, 4])
        True
        """
        tmp = self.list1[1:][:]
        tmp.sort(*kwargs)
        self.list1 = [0] + tmp