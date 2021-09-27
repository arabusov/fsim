from panda3d.core import loadPrcFileData

confVars = '''
window-title Forest Game
'''

loadPrcFileData('', confVars)

from direct.showbase.ShowBase import ShowBase
from panda3d.core import CollisionTraverser, CollisionHandlerPusher
from panda3d.core import CollisionNode, CollisionSphere, CollisionInvSphere


class MyGame(ShowBase):
    def __init__(self):
        ShowBase.__init__(self)

        self.disableMouse()
        
        self.cTrav = CollisionTraverser()
        
        pusher = CollisionHandlerPusher()
        pusher.setHorizontal(True)

        self.env = self.loader.loadModel("models/environment")
        self.env.setScale(0.25, 0.25, 0.25)
        self.env.setPos(-8, 42, -3)
        self.env.reparentTo(self.render)
        envBounds = self.env.getBounds()
        envCenter = envBounds.getCenter()
        envRad = envBounds.getRadius() * 0.7
        
        cNode = CollisionNode("environment")
        cNode.addSolid(CollisionInvSphere(envCenter, envRad))
        envC = self.render.attachNewNode(cNode)
        
        camBounds = self.camera.getBounds()
        camCenter = camBounds.getCenter()
        camRad = 5
         
        cNode = CollisionNode("camera")
        cNode.addSolid(CollisionSphere(camCenter, camRad))
        camC = self.camera.attachNewNode(cNode)
        camC.show()
        
        self.cTrav.addCollider(camC, pusher)
        pusher.addCollider(camC, self.camera)

        self.accept("w", self.__key, ["w"])
        self.accept("a", self.__key, ["a"])
        self.accept("s", self.__key, ["s"])
        self.accept("d", self.__key, ["d"])
        self.accept("arrow_left", self.__key, ["left"])
        self.accept("arrow_right", self.__key, ["right"])

    def __key(self, key):
        if key == "w":
            self.camera.setY(self.camera, 2)
        elif key == "a":
            self.camera.setX(self.camera, -2)
        elif key == "s":
            self.camera.setY(self.camera, -2)
        elif key == "d":
            self.camera.setX(self.camera, 2)
        elif key == "left":
            self.camera.setH(self.camera, 2)
        elif key == "right":
            self.camera.setH(self.camera, -2)


game = MyGame()
game.run()


