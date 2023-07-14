#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yzj
# @FileName :hhh.py
# @Time :2021/10/29 17:08

# 在python 中，列表变量用+=本质上是
# 在执行列表变量的extend 方法，不会修改变量的引用


class Gun:

    def __init__(self, model):
        self.model = model

        self.bullet_count = 0

    # def add_count(self):
    def add_count(self, count):
        self.bullet_count += count

    def shoot(self):
        # if self.bullet_count() <= 0:
        # 1.判断子弹数量
        if self.bullet_count <= 0:
            print("[%s] 还没有子弹..." % self.model)

            return

        # self.bullet_count() -= 1
        # 2.剩余子弹
        self.bullet_count -= 1

        # 3.发射子弹
        print("[%s] 突突突。。。[%s]" % (self.model, self.bullet_count))


class Soldier:

    def __init__(self, name):
        self.name = name

        self.gun = None

    def fire(self):
        # 1.判断士兵是否有枪
        if self.gun is None:
            print("[%s] 还没有抢。。。" % self.name)

            return

        # 2.高喊口号
        print("[%s] 冲啊。。。" % self.name)

        # 3.安装子弹
        self.gun.add_count(50)

        # 4.发射子弹
        self.gun.shoot()


# 1.创建枪对象
Ak47 = Gun("Ak47")

# 2.创建士兵对象
xusanduo = Soldier("许三多")

xusanduo.gun = Ak47
xusanduo.fire()
