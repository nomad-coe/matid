(self.webpackChunkdocs=self.webpackChunkdocs||[]).push([[20],{8019:(e,s,t)=>{e.exports={src:{srcSet:t.p+"assets/ideal-img/sbc.1c4fd99.640.png 640w,"+t.p+"assets/ideal-img/sbc.ff70971.1030.png 1030w",images:[{path:t.p+"assets/ideal-img/sbc.1c4fd99.640.png",width:640,height:398},{path:t.p+"assets/ideal-img/sbc.ff70971.1030.png",width:1030,height:640}],src:t.p+"assets/ideal-img/sbc.1c4fd99.640.png",toString:function(){return t.p+"assets/ideal-img/sbc.1c4fd99.640.png"},placeholder:void 0,width:640,height:398},preSrc:"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAGCAYAAAD68A/GAAAACXBIWXMAAC4jAAAuIwF4pT92AAAA4ElEQVR4nGWPP0tCUQBH7wdybOojOESfwKmhUWhoaQgSAm0OXIoKfOAfwsKX2Jg8FTNCCoxX0UuRRK/vVXJ9euOeqMHF33oOHH7CsDxjDE6tQS53Trv9wGgkEa7bwcocUbYvGI8HC7HrdXCqJbreC1rPEVe3NluHmyTPNnjr2aiJoi8lheYJmZs4radTgk+JaFXSHG9HyKZW8R4rfH2H+IHPq9fAuc/y/uGgfwziuVqgtLtGcSfKwL37T4eTgOuDGPn4CnVrH18OEX9A6xl6phZn5qGimU9ymVinZu0RThW/cM3OGKl7zVMAAAAASUVORK5CYII="}},865:(e,s,t)=>{"use strict";t.r(s),t.d(s,{assets:()=>m,contentTitle:()=>u,default:()=>p,frontMatter:()=>d,metadata:()=>h,toc:()=>g});var n=t(5893),r=t(1151),i=t(5944),a=t(8019),l=t.n(a),o=t(7416);const c="from matid.clustering import SBC\nimport ase.io\nfrom ase.visualize import view\n\nsystem = ase.io.read('system.xyz')\n\nsbc = SBC()\nclusters = sbc.get_clusters(system)\n\nfor cluster in clusters:\n    print(cluster.indices)\n    view(cluster.get_cell())\n",d={sidebar_position:3},u="Symmetry-based Clustering (SBC)",h={id:"learn/symmetry-based-clustering",title:"Symmetry-based Clustering (SBC)",description:"There are several tools for building atomistic systems from components, but the",source:"@site/docs/learn/symmetry-based-clustering.mdx",sourceDirName:"learn",slug:"/learn/symmetry-based-clustering",permalink:"/matid/docs/learn/symmetry-based-clustering",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:3,frontMatter:{sidebar_position:3},sidebar:"docsSidebar",previous:{title:"Installation",permalink:"/matid/docs/get-started/installation"},next:{title:"Symmetry Analysis",permalink:"/matid/docs/learn/symmetry-analysis"}},m={},g=[];function f(e){const s={a:"a",code:"code",h1:"h1",p:"p",...(0,r.a)(),...e.components};return(0,n.jsxs)(n.Fragment,{children:[(0,n.jsx)(s.h1,{id:"symmetry-based-clustering-sbc",children:"Symmetry-based Clustering (SBC)"}),"\n",(0,n.jsx)(s.p,{children:"There are several tools for building atomistic systems from components, but the\nreverse process of identifying these components from an existing system is much\nharder. MatID utilizes a custom clustering algorithm, called Symmetry-based\nClustering (SBC), which can cluster atoms based on local translational symmetry.\nEssentially this means that atoms which are built from the same underlying unit\ncell, will get clustered together."}),"\n",(0,n.jsx)(i.Z,{img:l(),style:{maxWidth:"35rem",marginRight:"auto",marginLeft:"auto"}}),"\n",(0,n.jsx)(s.p,{children:"Unlike clustering methods that work solely on the basis of atomic distances or\nchemical identities, SBC can produce meaningful clustering results even for\nstructures like grain boundaries and defected structures, where other tools\nwould struggle. SBC works for both finite and periodic structures, and can deal\nwith noise and curvature in the structures. SBC currently returns clusters where\nthe underlying unit cell is repeated either in two or three distinct directions\n(one-dimensional clusters, such as polymers are currently not handled)."}),"\n",(0,n.jsxs)(s.p,{children:["Clustering is performed using the ",(0,n.jsx)(s.a,{href:"/docs/reference/sbc",children:(0,n.jsx)(s.code,{children:"SBC"})})," class. The basic syntax is relatively\nsimple: you have to initialize the SBC class with parameters that are suitable\nfor your use case (sensible defaults are provided), and then call the\n",(0,n.jsx)(s.a,{href:"/docs/reference/sbc#get_clusters",children:(0,n.jsx)(s.code,{children:"get_clusters"})}),"-method. The following demonstrates this on a existing structure\nfile:"]}),"\n","\n",(0,n.jsx)(o.Z,{language:"python",children:`${c.split("\n").slice(0,8).join("\n")}`}),"\n",(0,n.jsxs)(s.p,{children:["The return value is a list of the found ",(0,n.jsx)(s.a,{href:"/docs/reference/cluster",children:(0,n.jsx)(s.code,{children:"Cluster"})})," instances. You can perform\nfurther analysis on these found clusters by using methods and attributes of this\nclass. E.g. to print out the indices of atoms belonging to the cluster and to\nvisualize the units cells of the found clusters, we could do the following:"]}),"\n",(0,n.jsx)(o.Z,{language:"python",children:`${c.split("\n").slice(9,12).join("\n")}`})]})}function p(e={}){const{wrapper:s}={...(0,r.a)(),...e.components};return s?(0,n.jsx)(s,{...e,children:(0,n.jsx)(f,{...e})}):f(e)}}}]);