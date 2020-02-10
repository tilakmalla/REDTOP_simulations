package REDTOPcontrib;
import OTPCClusterers.OTPCClumpFinder;
import hep.aida.IAnalysisFactory;
import hep.aida.ICloud1D;
import hep.aida.ICloud2D;
import hep.aida.ITree;
import hep.aida.ITuple;
import hep.aida.ITupleFactory;
import hep.physics.vec.BasicHep3Vector;
import org.lcsim.util.Driver;
import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import org.lcsim.conditions.ConditionsManager;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.geometry.Calorimeter;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.IDDecoder;
import org.lcsim.geometry.Subdetector;
import org.lcsim.geometry.util.DetectorLocator;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.hitmap.HitMap;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.event.LCIOParameters;
import org.lcsim.event.base.BaseTrack;
import org.lcsim.geometry.FieldMap;
import org.lcsim.recon.mcTrackFinder.MCTrack;
import org.lcsim.recon.tracking.seedtracker.TrackState;
import org.lcsim.spacegeom.SpacePoint;
import org.lcsim.spacegeom.SpaceVector;
import org.lcsim.util.swim.HelixSwimmer;
import org.lcsim.recon.cat.util.BasicTrack;
import org.lcsim.recon.tracking.trfutil.WAvg;

/**
 * *********************************************************************************
 * Driver to study the resolution of the OTPC in REDTOP. Requires that the hits
 * has been already digitized and the digi collections for the OTPCnur
 * (barrel+endcaps) be present in the event.
 
 * Steps in the drivers are as follow: Step 1 - general analysis of the particle
 * species present in the event. No ntuples are filled, since the info should be
 * used in other parts of the program.
 *
 * Step 2 - Study of the ring aperture. The hits (digi) found in the OTPC are
 * associated to the MCparticle which generated them. The angle wrt theredtop.200mev_90deg.v501.singlemuon1.opticalphotons.10001.1.00.slcio
 * direction of the particle (swum to the inner radius of the aerogel) is
 * calculated and compared to the expected Cerenkov angle. Two ntuples are
 * filled: on per each track in the event and one per each digi. Goal is to
 * compute the resolution on the Cerenkov angle
 *
 * Step 3 - Study of the impact point of a track in the sensors of the OTPC and
 * of the digi generated in the gas. The track is swum to the boundaries of the
 * OTPC. The impact point and the intersection of its direction with the OTPC
 * boundaries (center of the Cerenkov ring) are compared with the closest digi's
 
 * found. One ntuple is filled per track with the above info.
 * @modified by Malla 2018/01/20
 */
public class mallaOTPCResolutionDriver extends Driver {
    String version = "0.10";
    
    //AIDA histogram declarations
     
    private static AIDA aida = AIDA.defaultInstance();
    ICloud1D h1;

    private ITupleFactory tf;
    private ITuple tuple,tupleOp, tupleDigi, tupleEle, tupleDebug;
    private IAnalysisFactory af;
    private ITree tree;
    private String AidaTreeName;
    private ICloud1D c1dCount, digiCount, c1dCount1,c1dCount2;
    private ICloud2D c2dCount;

    protected int EventCount;
    protected int nPrinted, nPrintMax;
    Detector REDTOPdetector;
    double Bfield;
    double thisBfield;
    private boolean _bdebug = false;
    private Hep3Vector PrimaryVtx = new BasicHep3Vector(0., 0., 0.);
    
    private Detector REDTOPDetector;
    //a contructor 
    public mallaOTPCResolutionDriver() {
        ConditionsManager mgr = ConditionsManager.defaultInstance();
        REDTOPDetector = DetectorLocator.findDetector("REDTOP", 0, mgr);
        AidaTreeName = "OTPCResolution";
        /*
         * Initialize histograms
         */
        tree = AIDA.defaultInstance().tree();
        tree.cd("/");
        tree.mkdir(AidaTreeName);
        tree.cd(AidaTreeName);

        // nt-ple related initialization
        af = aida.analysisFactory();
        tf = af.createTupleFactory(aida.tree());
        String names = "float pe; float digi; float charge; float CerenkovAngle; float betaTrack; float momentumTrack; "
                + "float costhetaTrack; float #mothers; float CerAngleFromDigi; float CerAngleFromDigiPrimaryTrack; float CerAngleFromDigiResidual;"
                + "float phiTrack; float thetaFromDigiResidual; float phiFromDigiResidual; float #trackabledaughters; float pdg;"
                + "float xdigiAve; float ydigiAve; float zdigiAve; float avgthetaCerenkovfromDigi";
        
        tuple = tf.create("OTPCResolution_tracks", "properties of tracks", names, "");
        String names1 = "double cerenkovangle;" 
                        + "double modAo; double PDG; double modpatAo;modpm;"
                        +"double cerenkovangle1;double modAo1 ;double modpatAo1";
        tupleOp = tf.create("OTPCResolution_tracks1", "properties of tracks1", names1, "");
       

        String namesDigi = "float ThetaCerdelta; float CerenkovAngle; float thetaCerFromDigi; float pdg; float beta; "
                + "float momentumTrack; float thetaTrack;float Corrected_Energy; float #digi; float #mothers; float #trackabledaughters; "
                + "float xdigi; float ydigi; float zdigi; float thetadigi;"
                + "float xprojectedDigi;float yprojectedDigi;zprojectedDigi; float thetaCerenkovfromDigi";
                
        tupleDigi = tf.create("OTPCResolution_digi", "properties of OTPCdigi", namesDigi, "");

        String namesEle = "float TrackSipmdelta; float CerenkovAngle; float zDigiOfTrackAtSensors; float pdg; float betaTrack; "
                + "float momentumTrack; float costhetaTrack; float peOfTrackAtSensors; float #DigiOfTrackAtSensors; float #mothers;"
                + "float RingCenterDelta; float digiweptbycone; float transversemomentumTrack; "
                + "float trueTrackSweep; float trackSweepResolution";
        tupleEle = tf.create("OTPCResolution_ele", "properties of OTPC digi for electrons", namesEle, "");
        String namesDebug = "float TrackSipmdelta; float CerenkovAngle; float zHit; float pdg; float betaTrack; "
                + "float momentumTrack; float costhetaTrack; float pe; float clustersize; float #mothers;"
                + "float RingCenterDelta; float digisweptbycone; float transversemomentumTrack;"
                + "float trueTrackSweep; float trackSweepResolution";
        if (_bdebug) {
            tupleDebug = tf.create("OTPCResolution_debug", "properties of OTPC digi for electrons", namesDebug, "");
        }

        EventCount = 0;
        nPrinted = 0;
        nPrintMax = 100;
    }
    
    /**
     * Startup initialization
     */
    protected void startOfData() {
        System.out.println("OTPCResolutionDriver version " + version
                + " initialization -------------------");
        printConfig();
        /*
         * Initialize histograms
         */
        tree = AIDA.defaultInstance().tree();
        tree.cd("/");
      //tree.mkdir(AidaTreeName);
       //tree.cd(AidaTreeName);        
        AIDA.defaultInstance().tree().cd("/" + AidaTreeName);
        c1dCount = aida.cloud1D("Particle Count by Event");
        c1dCount1 = aida.cloud1D(" cerenkovanglep");
        c1dCount2 = aida.cloud1D("energyof p");
        digiCount = aida.cloud1D("OTPC digi Count by Event");
        c2dCount= aida.cloud2D("thetacerenkov vs radius");

        System.out.println("...end OTPCResolutionDriver initialization -------------------");
    }

    /**
     * Event loop
     *
     * @param event The event.
     */
    @Override
    protected void process(EventHeader event) {
        AIDA.defaultInstance().tree().cd("/" + AidaTreeName);

        EventCount++;
        REDTOPdetector = event.getDetector();
        FieldMap fmap = REDTOPdetector.getFieldMap();
        double[] origin = new double[]{0., 0., 0.};
        Bfield = fmap.getField(origin)[2];
        System.out.println("Bfield is:"+ Bfield);
      


    /* ***************************************************************************************************
         * First section of the analysis: look at the lists of particles present
         * in the event and update the counters.
     *****************************************************************************************************/
        int opticalPhotons = 0;
        int regularPhotons = 0;
        int regularCharged = 0;
        int etas = 0;
        int pi0s = 0;
        int picharged = 0;
        int muons = 0;
        int electrons = 0;
        int protons = 0;
        int neutrons = 0;
        int primarypi0 = 0;
        int primarypich = 0;
        int primaryelectrons = 0;
        int primarymuons = 0;
        int primaryphoton = 0;
        int primaryproton = 0;
        int primaryneutron = 0;
        int primarykch = 0;
        int primaryk0 = 0;
        int primaryHe3 = 0;

        int detectablepi0 = 0;
        int detectablepich = 0;
        int detectableelectrons = 0;
        int detectablemuons = 0;
        int detectablephoton = 0;
        int detectableproton = 0;
        int detectableneutron = 0;
        int detectablekch = 0;
        int detectablek0 = 0;
        int detectableHe3 = 0;

        // loop over MCParticles
        boolean passFilter = false;
        for (MCParticle p : event.getMCParticles()) {
            double energyOfp = p.getEnergy();
           List<MCParticle> opticalPhotonList = new ArrayList<MCParticle> (); 
            
            int gs = p.getGeneratorStatus();
            System.out.println(gs);
            MCParticle.SimulatorStatus ss = p.getSimulatorStatus();
            int pid = Math.abs(p.getPDGID());
            double pch = p.getCharge();
            boolean isopticalPhoton = (gs == 0 && (pid == 0) && p.getEnergy() < 1.0E-8);
            boolean isPrimary = (p.getParents().size() <= 0);
            boolean isTrackable = IsOTPCTrackable(p);
            boolean isShowerable = IsADRIANOShowerable(p);

            if (pid == 111 && isPrimary) // primary pi0
            {
                primarypi0++;
            }
            if (pid == 211 && isPrimary) // primary pi+
            {
                primarypich++;
            }
            if (pid == 11 && isPrimary) // primary e- 
            {
                primaryelectrons++;
            }
            if (pid == 13 && isPrimary) // primary mu+
            {
                primarymuons++;
            }
            if (pid == 22 && isPrimary) // primary gamma
            {
                primaryphoton++;
            }
            if (pid == 2212 && isPrimary) // primary proton
            {
                primaryproton++;
            }
            if (pid == 2112 && isPrimary) // primary neutron
            {
                primaryneutron++;
            }
            if (pid == 321 && isPrimary) // primary K+
            {
                primarykch++;
            }
            if (pid == 311 && isPrimary) // primary K0
            {
                primaryk0++;
            }
            if (pid == 1000020030 && isPrimary) // primary He3
            {
                primaryHe3++;
            }

            if (pid == 111 && isShowerable) // detectable pi0
            {
                detectablepi0++;
            }
            if (pid == 211 && isTrackable) // detectable pi+
            {
                detectablepich++;
            }
            if (pid == 11 && isTrackable) // detectable e- 
            {
                detectableelectrons++;
            }
            if (pid == 13 && isTrackable) // detectable mu+
            {
                detectablemuons++;
            }
            if (pid == 22 && isShowerable) // detectable gamma
            {
                detectablephoton++;
            }
            if (pid == 2212 && isShowerable) // detectable proton
            {
                detectableproton++;
            }
            if (pid == 2112 && isShowerable) // detectable neutron
            {
                detectableneutron++;
            }
            if (pid == 321 && isTrackable) // detectable K+
            {
                detectablekch++;
            }
            if (pid == 311 && isShowerable) // detectable K0
            {
                detectablek0++;
            }
            if (pid == 1000020030 && isShowerable) // detectable He3
            {
                detectableHe3++;
            }

            // secondary particle analysis
            if (!isopticalPhoton && gs == 0 && pid == 22) //  found a FS regular photon
            {
                regularPhotons++;
            }
            if (isopticalPhoton) //  found an optical photon
            {
                opticalPhotons++;
                opticalPhotonList.add(p);
            }
            boolean isPi0 = (gs == 0 && (pid == 111));
            if (isPi0) {
                pi0s++;   // this is a pi0
            }
            boolean isEta = (gs == 0 && (pid == 221));
            if (isEta) {
                etas++;   // this is a eta
            }
            boolean isNeutron = (gs == 0 && (pid == 2112));
            if (isNeutron) {
                neutrons++;   // this is a neutron
            }
            if (pch != 0.0 && gs == 0) //  found a stopped charged particles
            {
                // Found a stopped particle
                regularCharged++;
                if (pid == 211) {
                    picharged++;
                }
                if (pid == 13) {
                    muons++;
                }
                if (pid == 11) {
                    electrons++;
                }
                if (pid == 2212) {
                    protons++;
                }
            }
        }
        if (event.getMCParticles() == null) {
            c1dCount.fill(0);
            return;
        }
        c1dCount.fill(event.getMCParticles().size());
          
        /**
         * ****************************************************************************************
         * Second section of the analysis: look at the lists of digi from the
         * OTPC and associate a MC particle to it. A HashMap is filled
         * containing the association of the MCparticle and its digi.
     	*******************************************************************************************
         */
        List<SimCalorimeterHit> OTPCBarrelDigiList = event.get(SimCalorimeterHit.class, "OpticalTpcBarrelDigiHits");
        List<SimCalorimeterHit> OTPCEndcapDigiList = event.get(SimCalorimeterHit.class, "OpticalTpcEndcapDigiHits");
       
        List<SimCalorimeterHit> OTPCBarrelHitList = event.get(SimCalorimeterHit.class, "OpticalTpcBarrelHits");
        List<SimCalorimeterHit> OTPCEndcapHitList = event.get(SimCalorimeterHit.class, "OpticalTpcEndcapHits");
        
        List<SimCalorimeterHit> OTPCAllHitList = new ArrayList<SimCalorimeterHit>();
        OTPCAllHitList.addAll(OTPCBarrelHitList);
        OTPCAllHitList.addAll(OTPCEndcapHitList);
        
        List<SimCalorimeterHit> OTPCAllDigiList = new ArrayList<SimCalorimeterHit>();
        OTPCAllDigiList.addAll(OTPCBarrelDigiList);
        OTPCAllDigiList.addAll(OTPCEndcapDigiList);
        
        digiCount.fill(OTPCAllDigiList.size());
        Map<MCParticle, Set<Long>> usedHitsMap = new HashMap<MCParticle, Set<Long>>();
        
        // Make a map of OTPC hits
        Map<Long, SimCalorimeterHit> inputSimhitHitMapOTPC = new HashMap<Long, SimCalorimeterHit>();
        addSimHits(event, "OpticalTpcBarrelHits", inputSimhitHitMapOTPC);
        addSimHits(event, "OpticalTpcEndcapHits", inputSimhitHitMapOTPC);
//        Set<Long> IDset = inputSimhitHitMapOTPC.keySet();

        Map<Long, SimCalorimeterHit> inputSimhitMapOTPC_fromTracks = new HashMap<Long, SimCalorimeterHit>();
        Map<Long, SimCalorimeterHit> inputSimhitMapOTPCfromOpticalPhotons = new HashMap<Long, SimCalorimeterHit>();
        if (_bdebug) {
            // only hits with contributions from tracks (no opticalphotons). For debugging purposes.
            addSimHitsfromTracks(event, "OpticalTpcBarrelHits", inputSimhitMapOTPC_fromTracks);
            addSimHitsfromTracks(event, "OpticalTpcEndcapHits", inputSimhitMapOTPC_fromTracks);

            // only hits with contributions from tracks (only opticalphotons). For debugging pourposes.
            addSimHitsfromFromOpticalPhotons(event, "OpticalTpcBarrelHits", inputSimhitMapOTPCfromOpticalPhotons);
            addSimHitsfromFromOpticalPhotons(event, "OpticalTpcEndcapHits", inputSimhitMapOTPCfromOpticalPhotons);
        }

        // Make a map of OTPC digi
//        HitMap inputDigiMapOTPC = new HitMap(OTPCBarrelDigiList);
//        inputDigiMapOTPC.putAll(new HitMap(OTPCEndcapDigiList));
        Map<Long, SimCalorimeterHit> inputDigiMapOTPC = new HashMap<Long, SimCalorimeterHit>();
        addSimHits(event, "OpticalTpcBarrelDigiHits", inputDigiMapOTPC);
        addSimHits(event, "OpticalTpcEndcapDigiHits", inputDigiMapOTPC);

        
        
        /* The follwing codes implement the approximation that each optical photons are generated 
        at inner boundary of the aerogel.
         @author Malla
         */
        
        
         Loop over the all digi in the OTPC
        
        for(SimCalorimeterHit hit:OTPCAllHitList){
         int ncontribs = hit.getMCParticleCount();
            for (int i = 0; i < ncontribs; i++) {
                MCParticle mcp = hit.getMCParticle(i);
                int gs = mcp.getGeneratorStatus();
                int pid = Math.abs(mcp.getPDGID());
                boolean isopticalPhoton = (gs == 0 && (pid == 0) && mcp.getEnergy() < 1.0E-8);
                if(isopticalPhoton){
                    double sIn = 75.0;
                    double sOut = 110.0;
                    Hep3Vector Ao = mcp.getOrigin();
                    double modAo = Ao.magnitude();
                    Hep3Vector A = mcp.getEndPoint();
                    if(modAo<sOut&& modAo>sIn){
                        Hep3Vector AoA = VecOp.sub(mcp.getEndPoint(),mcp.getOrigin());
                        MCParticle opticalMother = mcp.getParents().get(0);
                    //double energymother = opticalMother.getEnergy();
                        Hep3Vector po = opticalMother.getOrigin();
                        Hep3Vector pm = opticalMother.getMomentum();
                        double modpm = pm.magnitude();
                    
                    
                    

                        BasicTrack basictrack = new BasicTrack();
                        basictrack.setHelixParameters(new SpacePoint(PrimaryVtx),new SpacePoint(po),new SpaceVector(pm), (int)opticalMother.getCharge());
                        HelixSwimmer swimmer = new HelixSwimmer(Bfield);
                        swimmer.setTrack(basictrack);


                        double tracklengthAtAo = swimmer.getTrackLengthToPoint(new SpacePoint(Ao));
                        SpaceVector patAo = swimmer.getMomentumAtLength(tracklengthAtAo);
                        double coscerenkovangle = VecOp.dot(AoA,patAo)/AoA.magnitude()/patAo.magnitude();
                        double cerenkovangle = Math.acos(coscerenkovangle);
                    
                        BaseTrack mcBaseTrack = new BaseTrack();
                        double[] pnt = {PrimaryVtx.x(),PrimaryVtx.y(),PrimaryVtx.z()};
                        mcBaseTrack.setReferencePoint(pnt);
                        mcBaseTrack.setRefPointIsDCA(false);
                        mcBaseTrack.setFitSuccess(false);
                        SpacePoint po1 = new SpacePoint(opticalMother.getOrigin());
                        SpaceVector momentumAtPo1 = new SpaceVector(opticalMother.getMomentum());
                        SpacePoint refPoint = new SpacePoint(PrimaryVtx);
                        LCIOParameters parameters = LCIOParameters.SpaceMomentum2Parameters(po1, momentumAtPo1, refPoint, (int)opticalMother.getCharge(), Bfield);
                        double[] fitPar = parameters.getValues();
                        mcBaseTrack.setTrackParameters(fitPar, Bfield);
                        HelixSwimmer swimmer1 = new HelixSwimmer(Bfield);
                        swimmer1.setTrack(mcBaseTrack);

                        double tracklengthAtAo1 = swimmer1.getTrackLengthToPoint(new SpacePoint(Ao));
                        SpaceVector patAo1 = swimmer1.getMomentumAtLength(tracklengthAtAo1);
                        double coscerenkovangle1 = VecOp.dot(AoA,patAo1)/AoA.magnitude()/patAo1.magnitude();
                        double cerenkovangle1 = Math.acos(coscerenkovangle1);
                    
                        tupleOp.fill(0, (double) cerenkovangle);
                        tupleOp.fill(1,(double)modAo);
                        tupleOp.fill(2,(double)opticalMother.getPDGID());
                        tupleOp.fill(3,(double)patAo.magnitude());
                        tupleOp.fill(4,(double)modpm);
                        tupleOp.fill(5,(double)cerenkovangle1);
                        tupleOp.fill(6,(double)tracklengthAtAo1);
                        tupleOp.fill(7,(double)patAo1.magnitude());
                        tupleOp.addRow();
                    }
                }
            }
        }        
        
           
            /**
             * ************************************************************************
             * Next piece of code is no longer needed since the Digi have been
             * recast to SimCalorimeterHits and the info related to contributing
             * MCparticle is available. In fact, that is the MCparticle which
             * started the chain ending into that digi. Code is retained for
             * reference or possible future use.
             *
             * SimCalorimeterHit shit; if
             * (inputSimhitHitMapOTPC.containsKey(cellID)) shit =
             * inputSimhitHitMapOTPC.get(cellID); else {
             * System.out.println("OTPCResolutionDriver: found digi with no
             * corresponding hit"); continue; } // --------> start loop on
             * contributions to hit int ncontribs = shit.getMCParticleCount();
             * for (int i=0; i<ncontribs; i++) { MCParticle mcp =
             * shit.getMCParticle(i); MCParticle mcpMother =
             * originalMotherMCP(mcp); Set<Long> usedHits =
             * usedHitsMap.get(mcpMother); if (usedHits == null) { // New used
             * hits list usedHits = new HashSet<Long>();
             * usedHitsMap.put(mcpMother, usedHits);
             *
             * } else { // Existing particle usedHits.add(hit.getCellID()); } }
             * // end loop on contributions
			**************************************************************************
             */
            // Save the hits in a HashMap
             
        
        for (SimCalorimeterHit hit : OTPCAllDigiList) {
            // Put this hit in the hitmap:
            Long cellID = hit.getCellID();
            int ncontribs = hit.getMCParticleCount();
            for (int i = 0; i < ncontribs; i++) {
                MCParticle mcp = hit.getMCParticle(i);
                //MCParticle mcpMother = originalMotherMCP(mcp);
                //Set<Long> usedHits = usedHitsMap.get(mcpMother);  // Save the hit and the primary particle which generated it
                Set<Long> usedHits = usedHitsMap.get(mcp);          // Save the hit and the  particle which generated it
                if (usedHits == null) {   // New used hits list
                    usedHits = new HashSet<Long>();
               //usedHitsMap.put(mcpMother, usedHits);
                    usedHitsMap.put(mcp, usedHits);
                    usedHits.add(hit.getCellID());
                } else {   // Existing particle
                    usedHits.add(hit.getCellID());
                }
            }   // end loop on contributions
        }   // end loop on OTPCAllDigiList

        
        
        /**
         * ****************************************************************************************
         * Third section of the analysis: look at the lists of MCparticles found
         * abov         e and add the logic for analyzing them and making ntuples. Ring
         * analysis: requires knowing where the tracks start to radiate. Aim is
         * to understand the uncertainty on the Cerenkov angle. Relevant for
         * muon and pion beta resolution.
     	 ******************************************************************************************
         */
     
         
        for(MCParticle mcp : usedHitsMap.keySet()) {   // loop over primary MC particles which made hits in the detector (and its daughter)
            
       
            
            MCParticle trackablemcp = recursiveTrackableMother(mcp);
            
            double beta = trackablemcp.getMomentum().magnitude() / trackablemcp.getEnergy();
            double CosCerenkovAngle = 1. / beta / nD;
            double CerenkovAngle = Math.acos(CosCerenkovAngle);
           
            Hep3Vector po = trackablemcp.getOrigin();
            Hep3Vector pm = trackablemcp.getMomentum();
            double magpm= pm.magnitude();
            BasicTrack basictrack = new BasicTrack();
            basictrack.setHelixParameters(new SpacePoint(PrimaryVtx),new SpacePoint(po),new SpaceVector(pm), (int)trackablemcp.getCharge());
            HelixSwimmer swimmer = new HelixSwimmer(Bfield);
            swimmer.setTrack(basictrack);
            
            double s = GetDistanceToInnerAerogel(swimmer);
            double sOut = GetDistanceToOuterAerogel(swimmer);
            double sAv= (s+sOut)/2;
            SpacePoint Ao = swimmer.getPointAtLength(sAv); 
            SpaceVector trackMomAtAerogel = swimmer.getMomentumAtLength(sAv);
            
            
            double costhetaTrack = VecOp.cosTheta(trackMomAtAerogel);
            double phiTrack = VecOp.phi(trackMomAtAerogel);
            int mothers = (trackablemcp.getParents().size() == 0 ? 0 : 1);
            


            BasicTrack superFastTrack = new BasicTrack();
            Hep3Vector superfasttrackMomAtAerogel = VecOp.mult(10000000000.0, trackMomAtAerogel);

            
            superFastTrack.setHelixParameters(new SpacePoint(PrimaryVtx),new SpacePoint(Ao),new SpaceVector(superfasttrackMomAtAerogel),-1);
            swimmer.setTrack(superFastTrack);
            double sFake = GetDistanceToOTPCSensors(swimmer);
            SpacePoint To = swimmer.getPointAtLength(sFake);
            Hep3Vector AoTo = VecOp.sub(To, Ao);


            // loop on digi to compute the measured Cerenkov angle and its resolution
            double pe = 0.0;
            WAvg avgthetaCerenkovfromDigi = new WAvg();
            WAvg avgCerAngle = new WAvg();
            WAvg avgCerAnglePrimaryTracksOnly = new WAvg();
            WAvg avgThetaPrimaryTracksOnly = new WAvg();
            WAvg avgPhiPrimaryTracksOnly = new WAvg();
            WAvg avgXPrimaryTracksOnly = new WAvg();
            WAvg avgYPrimaryTracksOnly = new WAvg();
            WAvg avgZPrimaryTracksOnly = new WAvg();
            Set<Long> cellIDSet = usedHitsMap.get(mcp);
            for (Long cellID : cellIDSet) {     // loop on digi associated to the trackable MCparticle
                SimCalorimeterHit digiOTPC = (_bdebug? inputSimhitMapOTPCfromOpticalPhotons.get(cellID) : inputDigiMapOTPC.get(cellID));
                if (digiOTPC == null) {
                    continue;
                }
                if (_bdebug) {
                    pe += digiOTPC.getRawEnergy();
                } else {
                    pe += digiOTPC.getCorrectedEnergy();
                }
                
                
     
     
                
                
/******************************************************* @author Malla***************************
                 * A plane is recovered by taking into account the equation of path of particle 
                 * hitting on the optical TPC and point of impact. Then the actual hit on 
                 * dodecagonal detector surface are projected on the place to get Cerenkov rings. 
**************************************************************************************************/
                
   
   
                // compute the vector between the point where the track hit the inner radius
                //of the aerogel and the sipm being struck by the Cerenkov photon
                
                Hep3Vector AoA = VecOp.sub(digiOTPC.getPositionVec(), Ao);
                double d = VecOp.dot(AoTo,AoTo)/VecOp.dot(AoA,AoTo);
             
                Hep3Vector vecAlongCerPhoton = VecOp.mult(d,AoA);
				
			
				// Position vector of point of hit projected on the plane
                Hep3Vector projectedDigi = VecOp.add(vecAlongCerPhoton,Ao);
                
                // Radius of cerenkov cone
                double radius = VecOp.sub(projectedDigi, To).magnitude();
                double base = AoTo.magnitude();
                double tanthetaCerfromDigi = radius/base;
                
               // semi-opening angle of cerenkov cone
                double thetaCerenkovfromDigi = Math.atan(tanthetaCerfromDigi);
                
                
                //original method
                double thetadigi = Math.acos(VecOp.cosTheta(AoA));
                double costhetadigi = VecOp.cosTheta(AoA);
                double phidigi = VecOp.phi(AoA);
                double costhetaCerFromDigi = VecOp.dot(AoA, trackMomAtAerogel) / AoA.magnitude() / trackMomAtAerogel.magnitude();
                double thetaCerFromDigi = Math.acos(costhetaCerFromDigi);
                double thetaCerDelta = thetaCerFromDigi - CerenkovAngle;
                double thetaCerDelta1 = thetaCerenkovfromDigi-CerenkovAngle;
                avgCerAngle.addPair(thetaCerFromDigi, digiOTPC.getCorrectedEnergy());
                avgthetaCerenkovfromDigi.addPair(thetaCerenkovfromDigi,digiOTPC.getCorrectedEnergy());
                
                
                // make an average of digi falling reasonably close to the ring
                double angletolerance = 0.070;  // 70 mrad
                if (mothers == 0 && Math.abs(thetaCerDelta) < angletolerance) {
                   
                    avgCerAnglePrimaryTracksOnly.addPair(thetaCerFromDigi, digiOTPC.getCorrectedEnergy());
                    avgThetaPrimaryTracksOnly.addPair(thetadigi, digiOTPC.getCorrectedEnergy());
                    avgPhiPrimaryTracksOnly.addPair(phidigi, digiOTPC.getCorrectedEnergy());
                    avgXPrimaryTracksOnly.addPair(digiOTPC.getPositionVec().x(), digiOTPC.getCorrectedEnergy());
                    avgYPrimaryTracksOnly.addPair(digiOTPC.getPositionVec().y(), digiOTPC.getCorrectedEnergy());
                    avgZPrimaryTracksOnly.addPair(digiOTPC.getPositionVec().z(), digiOTPC.getCorrectedEnergy());
                    
                    
                }
                
                               
               //fill ntuple
               
                tupleDigi.fill(0, (float) thetaCerDelta); 
                tupleDigi.fill(1, (float) CerenkovAngle);
                tupleDigi.fill(2, (float) thetaCerFromDigi);
                tupleDigi.fill(3, (float) Math.abs(trackablemcp.getPDGID()));
                tupleDigi.fill(4, (float) beta);
                tupleDigi.fill(5, (float) trackablemcp.getMomentum().magnitude());
                tupleDigi.fill(6, (float) (Math.acos(costhetaTrack)));
                tupleDigi.fill(7, (float) digiOTPC.getCorrectedEnergy());
                tupleDigi.fill(8, (float) cellIDSet.size());
                tupleDigi.fill(9, (float) mothers);
                tupleDigi.fill(11, (float) digiOTPC.getPositionVec().x());
                tupleDigi.fill(12, (float) digiOTPC.getPositionVec().y());
                tupleDigi.fill(13, (float) digiOTPC.getPositionVec().z());
                tupleDigi.fill(14, (float) thetadigi);
                tupleDigi.fill(15,(float) projectedDigi.x());
                tupleDigi.fill(16,(float) projectedDigi.y());
                tupleDigi.fill(17,(float) projectedDigi.z()); 
                tupleDigi.fill(18, (float) thetaCerenkovfromDigi);
                tupleDigi.addRow();
            }   // end loop on digi associated to the trackable MCparticle
            
            tuple.fill(0, (float) pe);
            tuple.fill(1, (float) cellIDSet.size());
            tuple.fill(2, (float) trackablemcp.getCharge());
            tuple.fill(3, (float) beta);
            tuple.fill(4, (float) CerenkovAngle);
            tuple.fill(5, (float) trackablemcp.getMomentum().magnitude());
            tuple.fill(6, (float) costhetaTrack);
            tuple.fill(7, (float) mothers);
            tuple.fill(8, (float) avgCerAngle.average());
            tuple.fill(9, (float) avgCerAnglePrimaryTracksOnly.average());
            tuple.fill(10, (float) (avgCerAnglePrimaryTracksOnly.average() - CerenkovAngle));
            tuple.fill(11, (float) phiTrack);
            tuple.fill(12, (float) (avgThetaPrimaryTracksOnly.average() - Math.acos(costhetaTrack)));
            tuple.fill(13, (float) (avgPhiPrimaryTracksOnly.average() - phiTrack));
            tuple.fill(15, (float) Math.abs(trackablemcp.getPDGID()));
            tuple.fill(16, (float) avgXPrimaryTracksOnly.average());
            tuple.fill(17, (float) avgYPrimaryTracksOnly.average());
            tuple.fill(18, (float) avgZPrimaryTracksOnly.average());
            tuple.fill(19,(float) avgthetaCerenkovfromDigi.average());
            tuple.addRow();
        } // end loop on tracks

       
       
       /**
         * ****************************************************************************************
         * Fourth section of the analysis: look at the lists of MCparticles
         * found above and add the logic for analyzing them and making ntuples.
         * Center of the ring analysis: requires knowing where the tracks start
         * to radiate and where they hit they the OTPC sensors. Aim is to
         * understand the uncertainty on the swipe of hits generated in the OTPC
         * gas. Relevant for electron momentum resolution.
     	 *******************************************************************************************
       */
        Subdetector detOTPCBarrel = REDTOPdetector.getSubdetector("OTPC_SIPMBarrel");  // OTPC_SIPMBarrel
        Calorimeter calOTPCBarrel = (Calorimeter) detOTPCBarrel;
        Subdetector detOTPCEndcap = REDTOPdetector.getSubdetector("OTPC_SIPMEndcaps");  // OTPC_SIPMEndcaps
        Calorimeter calOTPCEndcap = (Calorimeter) detOTPCEndcap;
        double sipmHole = calOTPCEndcap.getInnerRadius();

        for (MCParticle mcp : usedHitsMap.keySet()) {   // loop over primary MC particles which made hits in the detector (and its daughter)
            // Analize the first tackable daughter (that is the last particle which has its origin 
            // before the aerogel and creat a track out of it.
            MCParticle trackablemcp = recursiveTrackableDauther(mcp);
            MCParticle trackablemcpAncestor = recursiveTrackableMother(mcp);
            if (trackablemcp == null) {
                continue;
            }
            
            // Create a track and swimm it to the OTPC sensors
            MCTrack mcTrack = new MCTrack(trackablemcp);
            Hep3Vector po = trackablemcp.getOrigin();
            Hep3Vector pm = trackablemcp.getMomentum();
            HelixSwimmer swimmer = new HelixSwimmer(Bfield);
         
            swimmer.setTrack(mcTrack);
            double tracklengthToParticleOrigin = swimmer.getTrackLengthToPoint(new SpacePoint(po));
            SpaceVector newpAtParticleOrigin = swimmer.getMomentumAtLength(tracklengthToParticleOrigin);

            // Swim the track  to the OTPC sensors and look for hits nearby. These are the hits
            // directly produced by the track hitting the sensor plane.
            boolean isBarrel = isOTPCSensorsBarrel(swimmer);
            double s = GetDistanceToOTPCSensors(swimmer);
            
            if (s == -1000.0) {
                System.out.println("OTPCResolutionDriver: found mcp with pdg=" + trackablemcp.getPDGID()
                        + " wich cannot be swum to either the barrel nor the endcap of the OPC.");
                continue;
            }

            // Found something.
            Hep3Vector intercept = swimmer.getPointAtDistance(s);//?
            // Obtain the unit vector giving the tangent:
            // Need to swim a little further to reach the sensitive slice of the detector
            // (warning: this is sensor dependent and should be based on geometry rather than hard coded)
            SpacePoint trackPosAtSensors = swimmer.getPointAtLength(s);//?
            SpaceVector trackMomentumAtSensors = swimmer.getMomentumAtLength(s);
            Hep3Vector tangent = VecOp.unit(VecOp.sub(trackPosAtSensors, intercept));
            Hep3Vector momentumAtIntercept = VecOp.mult(pm.magnitude(), tangent);//??
			//TrackExtrapolationInfo info = new TrackExtrapolationInfo(intercept, momentumAtIntercept, trackablemcp, swimmer);

            TrackState aerogelState = new TrackState(trackPosAtSensors.x(), trackPosAtSensors.y(), trackPosAtSensors.z(), trackMomentumAtSensors.x(), trackMomentumAtSensors.y(), trackMomentumAtSensors.z());
            double costhetaTrack = VecOp.cosTheta(trackMomentumAtSensors);
            int mothers = (trackablemcp.getParents().size() == 0 ? 0 : 1);

            // Check if the particles ends up in the hole
            if (trackPosAtSensors.rxy() < sipmHole) {
                System.out.println("OTPCResolutionDriver: found mcp with pdg=" + mcp.getPDGID()
                        + " which ends up in the front hole at pos." + trackPosAtSensors.x() + " ; "
                        + trackPosAtSensors.y() + " ; " + trackPosAtSensors.z());
            }

            // Get the OTPC sensitive cell corresponding to where the track hit the sensors plane
            IDDecoder calOTPCdec = GetOTPCDecoder(calOTPCBarrel, calOTPCEndcap, trackPosAtSensors);
            if (calOTPCdec == null) {
                // the track is escaping the detector before hitting the ITPC sensors
                System.out.println("OTPCResolutionDriver: found mcp with pdg=" + mcp.getPDGID()
                        + " hitting no sipm.");
                continue;
            }
            long cellIDOTPC;
            try {
                cellIDOTPC = calOTPCdec.findCellContainingXYZ(trackPosAtSensors);
            } catch (Exception e) {
                System.out.println("OTPCResolutionDriver: ERROR - found mcp with pdg=" + mcp.getPDGID()
                        + " hitting a non sensitive module.");
                System.out.println(e.getMessage());
                continue;
            }
            long cellIDwithHit = 0;
            List<CalorimeterHit> neighbourTrackPosAtSensorsDigi = GetListOfNeighbourDigi(cellIDOTPC, calOTPCdec, inputDigiMapOTPC, 10);
            List<CalorimeterHit> neighbourTrackHitFromNonOpticalPhotons = new Vector<CalorimeterHit>();
            if (_bdebug) {
                neighbourTrackHitFromNonOpticalPhotons = GetListOfNeighbourDigi(cellIDOTPC, calOTPCdec, inputSimhitMapOTPC_fromTracks, 10);
            }

            if (neighbourTrackPosAtSensorsDigi.size() == 0) {
                System.out.println("OTPCResolutionDriver: found mcp with pdg=" + mcp.getPDGID()
                        + " with no OTPC digi near where the tracks intercept the OTPC sensors.");
                continue;
            } // cannot find an hit nearby where it is expected

            double closestDistance = 9999.0;
            int findex = -1;
            // loop over the found digi and select the one closest to the expected point
            for (int fi = 0; fi < neighbourTrackPosAtSensorsDigi.size(); fi++) {
                Hep3Vector digiV = neighbourTrackPosAtSensorsDigi.get(fi).getPositionVec();
                Hep3Vector diffV = VecOp.sub(digiV, trackPosAtSensors);
                if (diffV.magnitude() < closestDistance) {
                    closestDistance = diffV.magnitude();
                    findex = fi;
                }
            }   // search trough found digits
            CalorimeterHit closestDigiToTrackPosAtSensors = neighbourTrackPosAtSensorsDigi.get(findex);
            Hep3Vector closestDigiToTrackVector = closestDigiToTrackPosAtSensors.getPositionVec();
            // Next is the segment between the OTPC cell hit by the track and the closest digi found
            Hep3Vector diffVector = VecOp.sub(closestDigiToTrackVector, trackPosAtSensors);

            Hep3Vector diffVectorNO = new BasicHep3Vector();
            if (_bdebug && neighbourTrackHitFromNonOpticalPhotons.size() > 0) {
                double closestDistanceNO = 9999.0;
                int findexNO = -1;
                // loop over the found digi and select the one closest to the expected point
                for (int fiNO = 0; fiNO < neighbourTrackHitFromNonOpticalPhotons.size(); fiNO++) {
                    Hep3Vector digiVNO = neighbourTrackHitFromNonOpticalPhotons.get(fiNO).getPositionVec();
                    Hep3Vector diffVNO = VecOp.sub(digiVNO, trackPosAtSensors);
                    if (diffVNO.magnitude() < closestDistanceNO) {
                        closestDistanceNO = diffVNO.magnitude();
                        findexNO = fiNO;
                    }
                }   // search trough found digits
                CalorimeterHit closestDigiToTrackPosNO = neighbourTrackPosAtSensorsDigi.get(findexNO);
                Hep3Vector closestDigiToTrackVectorNO = closestDigiToTrackPosNO.getPositionVec();
                diffVectorNO = VecOp.sub(closestDigiToTrackVectorNO, trackPosAtSensors);
            }

            //Get the refractive index from compact.xml
            //Subdetector detAerogelOut = REDTOPdetector.getSubdetector("OTPC_AerogelOut");
            //Calorimeter calAerogelOut = (Calorimeter)detAerogelOut;
            //double innerR = calAerogelOut.getDetectorElement().getGeometry()
            double nD = 1.000444;   // refractive index (hard coded for now) 
            double beta = trackablemcp.getMomentum().magnitude() / trackablemcp.getEnergy();
            double CosCerenkovAngle = 1. / beta / nD;
            double CerenkovAngle = (CosCerenkovAngle <= 1.0 ? Math.acos(CosCerenkovAngle) : -1.0);    // for particles other than electrons

            // Now, repeat the analysis to search for the center of the ring generated
            // by the aerogel. Thiis is the intercept between the direction of the track at the aerogel
            // and the OTPC sensors. Do that by creating a fake neutral track which will be swum to
            // the OTPC sensors.
            // First, swim the current track to the Aerogel
            double sAG = GetDistanceToInnerAerogel(swimmer);
            SpacePoint trackPosAtAerogel = swimmer.getPointAtLength(sAG);
            SpaceVector trackMomAtAerogel = swimmer.getMomentumAtLength(sAG);
            // Second, create a fake superfast track (wich approximate a straight line, 
            // and swim it to the OTPCsensors
            BasicTrack superfastfaketrack = new BasicTrack();
            Hep3Vector superfasttrackMomAtAerogel = VecOp.mult(1000000000000000.0, trackMomAtAerogel);
            superfastfaketrack.setHelixParameters(new SpacePoint(new BasicHep3Vector()), new SpacePoint(trackPosAtAerogel),
                    new SpaceVector(superfasttrackMomAtAerogel), 1);
            swimmer.setTrack(superfastfaketrack);
            double sFake = GetDistanceToOTPCSensors(swimmer);
            if (sFake == -1000.0) {
                System.out.println("OTPCResolutionDriver: found fake neutral trak with pdg=" + trackablemcp.getPDGID()
                        + " wich cannot be swum to either the barrel nor the endcap of the OPC.");
                continue;
            }

            // Found something.
            SpacePoint faketrackPosAtSensors = swimmer.getPointAtLength(sFake);
            SpaceVector faketrackMomentumAtSensors = swimmer.getMomentumAtLength(sFake);
            IDDecoder RingCenterdec = GetOTPCDecoder(calOTPCBarrel, calOTPCEndcap, faketrackPosAtSensors);
            if (calOTPCdec == null) {
                // the track is escaping the detector before hitting the ITPC sensors
                System.out.println("OTPCResolutionDriver: found mcp with pdg=" + mcp.getPDGID()
                        + " hitting no sipm.");
                continue;
            }
            long cellIDRingCenter;
            try {
                cellIDRingCenter = RingCenterdec.findCellContainingXYZ(faketrackPosAtSensors);
            } catch (Exception e) {
                System.out.println("OTPCResolutionDriver: ERROR - found mcp with pdg=" + mcp.getPDGID()
                        + " with direction not intercepting an OTPC sensitive module.");
                System.out.println(e.getMessage());
                continue;
            }
            long cellIDRingCenterwithHit = 0;
            List<CalorimeterHit> neighbourRingCenterDigi = GetListOfNeighbourDigi(cellIDRingCenter, RingCenterdec, inputDigiMapOTPC, 10);
            if (neighbourRingCenterDigi.size() == 0) {
                System.out.println("OTPCResolutionDriver: found mcp with pdg=" + mcp.getPDGID()
                        + " with no OTPC digi near the center of the ring.");
               continue;
            } // cannot find an hit nearby where it is expected

            double closestDistanceToRingCenter = 9999.0;
            int cindex = -1;
            // loop over the found digi and select the one closest to the expected point
            for (int fi = 0; fi < neighbourRingCenterDigi.size(); fi++) {
                Hep3Vector digiV = neighbourRingCenterDigi.get(fi).getPositionVec();
                // Next is the segment between the center of the ring and the closest OTPC digi
                Hep3Vector diffV = VecOp.sub(digiV, faketrackPosAtSensors);
                if (diffV.magnitude() < closestDistanceToRingCenter) {
                    closestDistanceToRingCenter = diffV.magnitude();
                    cindex = fi;
                }
            }   // search trough found digits
            double SipmSegmentbeweenRingCenterandTrackHit = -100.0;
            double SegmentbeweenRingCenterandTrackPosAtsensor = -100.0;
            if (cindex > -1) {
                CalorimeterHit closestDigiToRingCenter = neighbourRingCenterDigi.get(cindex);
                Hep3Vector closestDigiToRingCenterVector = closestDigiToRingCenter.getPositionVec();
                Hep3Vector diffVectorRingCenter = VecOp.sub(closestDigiToRingCenterVector, faketrackPosAtSensors);
                Hep3Vector DigiSweep = VecOp.sub(closestDigiToRingCenterVector, closestDigiToTrackVector);
                SipmSegmentbeweenRingCenterandTrackHit = DigiSweep.magnitude();
                SegmentbeweenRingCenterandTrackPosAtsensor = diffVectorRingCenter.x();
            }
            // Next is the segment between the true ring center and the true OTPC point struck by the track
            Hep3Vector TrueSweep = VecOp.sub(faketrackPosAtSensors, trackPosAtSensors);

            // Next is the segment between the digi closest to the ring center and the digi closest 
            // to the sipm struck by the track
           // Hep3Vector DigiSweep = new BasicHep3Vector(0.,0.,0.);
           //if (cindex>-1) 
           // DigiSweep = VecOp.sub(closestDigiToRingCenterVector, closestDigiToTrackVector);
            double trackPt = Math.sqrt(trackablemcp.getPX() * trackablemcp.getPX()
                    + trackablemcp.getPY() * trackablemcp.getPY());
            
			// Fill the ntuples
            tupleEle.fill(0, (float) diffVector.x());
            tupleEle.fill(1, (float) CerenkovAngle);
            tupleEle.fill(2, (float) closestDigiToTrackPosAtSensors.getPositionVec().z());
            tupleEle.fill(3, (float) Math.abs(trackablemcp.getPDGID()));
            tupleEle.fill(4, (float) beta);
            tupleEle.fill(5, (float) trackablemcp.getMomentum().magnitude());
            tupleEle.fill(6, (float) costhetaTrack);
            tupleEle.fill(7, (float) closestDigiToTrackPosAtSensors.getCorrectedEnergy());
            tupleEle.fill(8, (float) neighbourTrackPosAtSensorsDigi.size());
            tupleEle.fill(9, (float) mothers);
            tupleEle.fill(10, (float) SegmentbeweenRingCenterandTrackPosAtsensor);
            tupleEle.fill(11, (float) SipmSegmentbeweenRingCenterandTrackHit);
            tupleEle.fill(12, (float) trackPt);
            tupleEle.fill(13, (float) TrueSweep.magnitude());
            tupleEle.fill(14, (float) (TrueSweep.magnitude() - SipmSegmentbeweenRingCenterandTrackHit));
            tupleEle.addRow();
            if (_bdebug) {
                // Fill the debug ntuples
                tupleDebug.fill(0, (float) diffVectorNO.x());
                tupleDebug.fill(1, (float) CerenkovAngle);
                tupleDebug.fill(2, (float) closestDigiToTrackPosAtSensors.getPositionVec().z());
                tupleDebug.fill(3, (float) Math.abs(trackablemcp.getPDGID()));
                tupleDebug.fill(4, (float) beta);
                tupleDebug.fill(5, (float) trackablemcp.getMomentum().magnitude());
                tupleDebug.fill(6, (float) costhetaTrack);
                tupleDebug.fill(7, (float) closestDigiToTrackPosAtSensors.getCorrectedEnergy());
                tupleDebug.fill(8, (float) neighbourTrackPosAtSensorsDigi.size());
                tupleDebug.fill(9, (float) mothers);
                tupleDebug.fill(10, (float) SegmentbeweenRingCenterandTrackPosAtsensor);
                tupleDebug.fill(11, (float) SipmSegmentbeweenRingCenterandTrackHit);
                tupleDebug.fill(12, (float) trackPt);
                tupleDebug.fill(13, (float) TrueSweep.magnitude());
                tupleDebug.fill(14, (float) (TrueSweep.magnitude() - SipmSegmentbeweenRingCenterandTrackHit));
                tupleDebug.addRow();
            }
        } // end loop on tracks


////////////////////////////////////////////////////////////////////////////////////////////
        if (nPrinted == nPrintMax) {
            nPrinted++;
            System.out.println("EventTrigger printout terminated: max print count"
                    + " (" + nPrintMax + " lines) reached");
        }

    }

    /**
     * Called by the framework when event processing is suspended, usually at
     * the end of data processing
     */
    protected void suspend() {
        System.out.println("OTPCResolutionDriver version " + version
                + " summary -----------------------------");
        printConfig();
        System.out.println("...EventCount = " + EventCount);
        System.out.println("...end OTPCResolutionDriver summary -----------------------------");
    }

    protected void printConfig() {
        System.out.println("...nPrintMax = " + nPrintMax + " lines");
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected boolean IsOTPCTrackable(MCParticle mcp) {
        double startZ = mcp.getOriginZ();
        double startR = mcp.getOriginX() * mcp.getOriginX() + mcp.getOriginY() * mcp.getOriginY();
        startR = sqrt(startR);

        Hep3Vector iv = mcp.getOrigin();
        Hep3Vector rv = mcp.getEndPoint();
        double rmag = rv.magnitude();
        double endZ = rv.z();
        double endR = rv.x() * rv.x() + rv.y() * rv.y();
        endR = sqrt(endR);

        Subdetector detAerogelIn = REDTOPdetector.getSubdetector("OTPC_AerogelIn");  //BeamPipe
        Calorimeter calAerogelIn = (Calorimeter) detAerogelIn;
        double innerR = calAerogelIn.getInnerRadius();
//        double innerR = calBeamPipeGas.getInnerRadius();
//        double outerR = calBeamPipeGas.getOuterRadius();            
//        double halfZ = calBeamPipeGas.getZLength() / 2;
//        double zlength = calBeamPipeGas.getZLength();
//        double  global_det_z1 = calBeamPipeGas.getInnerZ();
//        double  global_det_z2 = calBeamPipeGas.getOuterZ();
//        boolean isBeamPipeGas = detBeamPipeGas.getDetectorElement().getGeometry().inside(iv) == Inside.INSIDE;
        Subdetector detAerogelOut = REDTOPdetector.getSubdetector("OTPC_AerogelOut");
        Calorimeter calAerogelOut = (Calorimeter) detAerogelOut;
        double outerR = calAerogelOut.getOuterRadius();
        double halfZ = calAerogelOut.getZLength() / 2;
        double zlength = calAerogelOut.getZLength();
        double global_det_z1 = calAerogelOut.getInnerZ();
        double global_det_z2 = calAerogelOut.getOuterZ();
//        boolean isAerogelOut = detBeamPipeGas.getDetectorElement().getGeometry().inside(iv) == Inside.INSIDE;

        if (startR < outerR && Math.abs(startZ) < halfZ
                && endR >= innerR) {
            return true;
        }

//        if (isAerogelOut) return true;
        return false;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected boolean IsADRIANOShowerable(MCParticle mcp) {
        double startZ = mcp.getOriginZ();
        double startR = mcp.getOriginX() * mcp.getOriginX() + mcp.getOriginY() * mcp.getOriginY();
        startR = sqrt(startR);

        Hep3Vector iv = mcp.getOrigin();
        Hep3Vector rv = mcp.getEndPoint();
        double rmag = rv.magnitude();
        double endZ = rv.z();
        double endR = rv.x() * rv.x() + rv.y() * rv.y();
        endR = sqrt(endR);

        Subdetector detADRIANOBarrel = REDTOPdetector.getSubdetector("ADRIANO_Barrel");  //ADRIANO barrel
        Calorimeter calADRIANOBarrel = (Calorimeter) detADRIANOBarrel;
        double innerR = calADRIANOBarrel.getInnerRadius();
        double outerR = calADRIANOBarrel.getOuterRadius();
        double halfZ = calADRIANOBarrel.getZLength() / 2;
        double zlength = calADRIANOBarrel.getZLength();
        double global_det_z1 = calADRIANOBarrel.getInnerZ();
        double global_det_z2 = calADRIANOBarrel.getOuterZ();

        if (startR < innerR && Math.abs(startZ) < halfZ) {
            return true;
        }
        return false;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected List<Cluster> findOTPCClumps(HitMap unusedHits) {
        //Clusterer clumpfinder = new OTPCClumpFinder();
        OTPCClumpFinder clumpfinder = new OTPCClumpFinder();
        List<Cluster> clumpClusters = clumpfinder.createClusters(unusedHits);
        for (Cluster clump : clumpClusters) {
            if (clump.getCalorimeterHits().size() == 0) {
                throw new AssertionError("clump has no hits");
            }
            if (clump.getCalorimeterHits().contains(null)) {
                throw new AssertionError("null hit in clump");
            }
            //   for (CalorimeterHit hit : clump.getCalorimeterHits()) {
            //       hitsInTree.remove(hit.getCellID());
            //   }
        }
        return clumpClusters;
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////
    private void addSimHits(EventHeader event, String name, Map<Long, SimCalorimeterHit> hitMap) {
        if (event.hasCollection(SimCalorimeterHit.class, name)) {
            List<SimCalorimeterHit> hitsToAdd = event.get(SimCalorimeterHit.class, name);
            for (SimCalorimeterHit hit : hitsToAdd) {
                hitMap.put(hit.getCellID(), hit);
            }
        }
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    private void addSimHitsfromTracks(EventHeader event, String name, Map<Long, SimCalorimeterHit> hitMap) {
        if (event.hasCollection(SimCalorimeterHit.class, name)) {
            List<SimCalorimeterHit> hitsToAdd = event.get(SimCalorimeterHit.class, name);
            for (SimCalorimeterHit hit : hitsToAdd) {
                // --------> start loop on contributions to hit
                boolean isFromTrack = false;
                for (int i = 0; i < hit.getMCParticleCount(); i++) {
                    MCParticle mcp = hit.getMCParticle(i);
                    //if (isFromEtaorPi0(mcp)) contribfrompi0oreta++;
                    double ei = hit.getContributedEnergy(i);
                    if (ei > 1.0e-8) {
                        isFromTrack = true;
                        break;
                    }
                }
                if (isFromTrack) {
                    hitMap.put(hit.getCellID(), hit);
                }
            }
        }
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    private void addSimHitsfromFromOpticalPhotons(EventHeader event, String name, Map<Long, SimCalorimeterHit> hitMap) {
        if (event.hasCollection(SimCalorimeterHit.class, name)) {
            List<SimCalorimeterHit> hitsToAdd = event.get(SimCalorimeterHit.class, name);
            for (SimCalorimeterHit hit : hitsToAdd) {
                // --------> start loop on contributions to hit
                boolean isFromOpticalPhoton = false;
                for (int i = 0; i < hit.getMCParticleCount(); i++) {
                    MCParticle mcp = hit.getMCParticle(i);
                    //if (isFromEtaorPi0(mcp)) contribfrompi0oreta++;
                    double ei = hit.getContributedEnergy(i);
                    if (ei < 1.0e-8) {
                        isFromOpticalPhoton = true;
                        break;
                    }
                }
                if (isFromOpticalPhoton) {
                    hitMap.put(hit.getCellID(), hit);
                }
            }
        }
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected MCParticle originalMotherMCP(MCParticle theParticle) {
        // look into up to five generations to find the primary mother
        MCParticle theOriginatorMCP = theParticle;

        // check if it is already a primary
        List<MCParticle> theGrandMothermcl = theParticle.getParents();
        if (theGrandMothermcl.size() == 0) {
            return theOriginatorMCP;
        }

        // found the mother #1
        MCParticle theGrandMother = theParticle.getParents().get(0);
        //theOriginatorPDG = theGrandMother.getPDGID();
        List<MCParticle> theGrandGrandMothermcl = theGrandMother.getParents();
        if (theGrandGrandMothermcl.size() == 0) {
            return theGrandMother;
        }

        // found the mother #2
        MCParticle theGrandGrandMother = theGrandMother.getParents().get(0);
        //theOriginatorPDG = theGrandGrandMother.getPDGID();
        List<MCParticle> theGrandGrandGrandMothermcl = theGrandGrandMother.getParents();
        if (theGrandGrandGrandMothermcl.size() == 0) {
            return theGrandGrandMother;
        }

        // found the mother #3
        MCParticle theGrandGrandGrandMother = theGrandGrandMother.getParents().get(0);
        //theOriginatorPDG = theGrandGrandGrandMother.getPDGID();
        List<MCParticle> theGrandGrandGrandGrandMothermcl = theGrandGrandGrandMother.getParents();
        if (theGrandGrandGrandGrandMothermcl.size() == 0) {
            return theGrandGrandGrandMother;
        }

        // found the mother #4
        MCParticle theGrandGrandGrandGrandMother = theGrandGrandGrandMother.getParents().get(0);
        //theOriginatorPDG = theGrandGrandGrandGrandMother.getPDGID();
        List<MCParticle> theGrandGrandGrandGrandGrandMothermcl = theGrandGrandGrandGrandMother.getParents();
        if (theGrandGrandGrandGrandGrandMothermcl.size() == 0) {
            return theGrandGrandGrandGrandMother;
        }

        // found no mother
        return null;
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected MCParticle recursiveTrackableDauther(MCParticle theParticle) {
        // look recursively  to look for a trackabel particle
        MCParticle theOriginatorMCP = theParticle;
        if (IsOTPCTrackable(theParticle)&&(theParticle.getCharge()!=0.0)) {
            return theOriginatorMCP;
        }

        // check if it is already a primary
        List<MCParticle> theDaugthers = theParticle.getDaughters();
        if (theDaugthers.size() == 0) {
            return null;     // found no trackable daughter
        }
        // Not trackable: recursively find the first trackable daughter
        for (MCParticle part : theDaugthers) {
            return recursiveTrackableDauther(part);
        }

        return null;

    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected MCParticle recursiveTrackableMother(MCParticle theParticle) {
        // look recursively  to look for a trackabel particle
        MCParticle theOriginatorMCP = theParticle;
        if (IsOTPCTrackable(theParticle)&& (theParticle.getCharge()!=0.0)) {
            return theOriginatorMCP;
        }

        // check if it is already a primary
        List<MCParticle> theParents = theParticle.getParents();
        if (theParents.size() == 0) {
            return null;     // found no trackable mother
        }
        for (MCParticle part : theParents) {
            return recursiveTrackableMother(part);
        }

        return null;

    }
    
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected boolean isOTPCSensorsBarrel(HelixSwimmer swimmer) {
        Subdetector detOTPCBarrel = REDTOPdetector.getSubdetector("OTPC_SIPMBarrel");  // OTPC_SIPMBarrel
        Calorimeter calOTPCBarrel = (Calorimeter) detOTPCBarrel;
        //double innerR = calBeamPipeGas.getInnerRadius();
        //double outerR = calOTPCAerogelIn.getOuterRadius(); 
        double sipmRadius = calOTPCBarrel.getInnerRadius();
        double delta0 = calOTPCBarrel.getDistanceToSensor(0);
        double deltaBarrel = detOTPCBarrel.getLayering().getDistanceToLayerSensorMid(0);
        int sipmNumSides = calOTPCBarrel.getNumberOfSides();
        Subdetector detOTPCEndcap = REDTOPdetector.getSubdetector("OTPC_SIPMEndcaps");  // OTPC_SIPMEndcaps
        Calorimeter calOTPCEndcap = (Calorimeter) detOTPCEndcap;
        //double sipmZ = calOTPCEndcap.getInnerZ();
        double sipmHole = calOTPCEndcap.getInnerRadius();
        double sipmZ = detOTPCEndcap.getLayering().getDistanceToLayerSensorMid(0);
        double deltaEndcap = 0.0;
        double s;
        double sZ = swimmer.getDistanceToZ(sipmZ);
        double sR = swimmer.getDistanceToPolyhedra(sipmRadius, sipmNumSides);
        boolean isBarrel = false;
        if (Double.isNaN(sR)) {   // track hits the endcap
            s = sZ;
            isBarrel = false;
        } else if (Double.isNaN(sZ)) {
            s = sR;
            isBarrel = true;
        } else if (swimmer.getDistanceToZ(sipmZ) < swimmer.getDistanceToPolyhedra(sipmRadius, sipmNumSides)) {   // endcap is closer
            s = sZ;
            isBarrel = false;
        } else if (swimmer.getDistanceToZ(sipmZ) >= swimmer.getDistanceToPolyhedra(sipmRadius, sipmNumSides)) {   // barrel is closer
            s = sR;
            isBarrel = true;
        } else {   // somethingh is wrong
            s = -1000.0;
            return isBarrel;
        }
        double delta = (isBarrel ? deltaBarrel : deltaEndcap);
        if (s < 0) {
            delta *= -1.0;
        }

        return isBarrel;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected double GetDistanceToOTPCSensors(HelixSwimmer swimmer) {
        Subdetector detOTPCBarrel = REDTOPdetector.getSubdetector("OTPC_SIPMBarrel");  // OTPC_SIPMBarrel
      
        Calorimeter calOTPCBarrel = (Calorimeter) detOTPCBarrel;
       //double innerR = calBeamPipeGas.getInnerRadius();
        //double outerR = calOTPCAerogelIn.getOuterRadius(); 
        double sipmRadius = calOTPCBarrel.getInnerRadius();
       double delta0 = calOTPCBarrel.getDistanceToSensor(0);
        double deltaBarrel = detOTPCBarrel.getLayering().getDistanceToLayerSensorMid(0);
        int sipmNumSides = calOTPCBarrel.getNumberOfSides();
        Subdetector detOTPCEndcap = REDTOPdetector.getSubdetector("OTPC_SIPMEndcaps");  // OTPC_SIPMEndcaps
        Calorimeter calOTPCEndcap = (Calorimeter) detOTPCEndcap;
        //double sipmZ = calOTPCEndcap.getInnerZ();
        double sipmHole = calOTPCEndcap.getInnerRadius();
        double sipmZ = detOTPCEndcap.getLayering().getDistanceToLayerSensorMid(0);
        double deltaEndcap = 0.0;
        double s;
        double sZ = swimmer.getDistanceToZ(sipmZ);
        double sR = swimmer.getDistanceToPolyhedra(sipmRadius, sipmNumSides);
        boolean isBarrel = false;
        if (Double.isNaN(sR)) {   // track hits the endcap
            s = sZ;
            isBarrel = false;
        } else if (Double.isNaN(sZ)) {
            s = sR;
            isBarrel = true;
        } else if (swimmer.getDistanceToZ(sipmZ) < swimmer.getDistanceToPolyhedra(sipmRadius, sipmNumSides)) {   // endcap is closer
            s = sZ;
            isBarrel = false;
        } else if (swimmer.getDistanceToZ(sipmZ) >= swimmer.getDistanceToPolyhedra(sipmRadius, sipmNumSides)) {   // barrel is closer
            s = sR;
            isBarrel = true;
        } else {   // somethingh is wrong
            s = -1000.0;
            return s;
        }
        double delta = (isBarrel ? deltaBarrel : deltaEndcap);
        if (s < 0) {
            delta *= -1.0;
        }

        return s + delta;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected double GetDistanceToInnerAerogel(HelixSwimmer swimmer) {
        // Get the inner aerogel for swimming tracks to it.
        Subdetector detOTPCAerogelIn = REDTOPdetector.getSubdetector("OTPC_AerogelIn");  // AerogelIn
        Calorimeter calOTPCAerogelIn = (Calorimeter) detOTPCAerogelIn;
       //double innerR = calBeamPipeGas.getInnerRadius();
       //double outerR = calOTPCAerogelIn.getOuterRadius(); 
        double otpcRadius = calOTPCAerogelIn.getInnerRadius();
        int otpcNumSides = calOTPCAerogelIn.getNumberOfSides();
        double otpcZ = calOTPCAerogelIn.getInnerZ();
        // Get the inner aerogel for swimming tracks to it.
        double s;
        double sZ = swimmer.getDistanceToZ(otpcZ);
        double sR = swimmer.getDistanceToPolyhedra(otpcRadius, otpcNumSides);
        if (Double.isNaN(sR)) {
            s = sZ;
        } else if (Double.isNaN(sZ)) {
            s = sR;
        } else {
            s = Math.min(swimmer.getDistanceToZ(otpcZ), swimmer.getDistanceToPolyhedra(otpcRadius, otpcNumSides));
        }
        return s;
    }
     protected double GetDistanceToOuterAerogel(HelixSwimmer swimmer) {
        // Get the inner aerogel for swimming tracks to it.
        Subdetector detOTPCAerogelOut = REDTOPdetector.getSubdetector("OTPC_AerogelOut");  // AerogelIn
        Calorimeter calOTPCAerogelOut = (Calorimeter) detOTPCAerogelOut;
        //double innerR = calBeamPipeGas.getInnerRadius();
        //double outerR = calOTPCAerogelIn.getOuterRadius(); 
        double otpcRadius = calOTPCAerogelOut.getOuterRadius();
        int otpcNumSides = calOTPCAerogelOut.getNumberOfSides();
        double otpcZ = calOTPCAerogelOut.getInnerZ();
        // Get the inner aerogel for swimming tracks to it.
        double s;
        double sZ = swimmer.getDistanceToZ(otpcZ);
        double sR = swimmer.getDistanceToPolyhedra(otpcRadius, otpcNumSides);
        if (Double.isNaN(sR)) {
            s = sZ;
        } else if (Double.isNaN(sZ)) {
            s = sR;
        } else {
            s = Math.min(swimmer.getDistanceToZ(otpcZ), swimmer.getDistanceToPolyhedra(otpcRadius, otpcNumSides));
        }
        return s;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected List<CalorimeterHit> GetListOfNeighbourDigi(long cellIDOTPC, IDDecoder calOTPCdec,
            Map<Long, SimCalorimeterHit> inputDigiMapOTPC, int deltaThetaPhi) {
        List<CalorimeterHit> neighbourTrackDigi = new Vector<CalorimeterHit>();

        //now, look for possible hits nearby: this is the detected position where the track hit the OTPC sensors
        calOTPCdec.setID(cellIDOTPC);
        if (!calOTPCdec.supportsNeighbours()) {
            throw new AssertionError("Can't get neighbours!");
        } else {
            // Parameters (hard-coded for now);
            int deltaLayer = 0;
            int deltaTheta = deltaThetaPhi;
            int deltaPhi = deltaThetaPhi;
            // Find all hits from [unusedHits] within range:
            long[] neighbourIDs = calOTPCdec.getNeighbourIDs(deltaLayer, deltaTheta, deltaPhi);
            // Look for cells which
            //   (a) are neighbours
            //   (b) are in the list "unusedHit", i.e. have a hit and aren't already in a clump/track/etc
            //   (c) aren't in the set "hitsUsedInClump", to avoid double-counting
            for (long neighbourID : neighbourIDs) { // loop over neighbours

                CalorimeterHit foundNeighbourDigi = inputDigiMapOTPC.get(neighbourID);// search OTPC digits
                if (foundNeighbourDigi != null) {
                           // cellIDwithHit = neighbourID;
                    neighbourTrackDigi.add(foundNeighbourDigi);
                }
            }   // end loop over neighbours
        }   // end if support neighbours
        return neighbourTrackDigi;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////

    protected IDDecoder GetOTPCDecoder(Calorimeter calOTPCBarrel, Calorimeter calOTPCEndcap, Hep3Vector trackPosAtSensors) {
        IDDecoder calOTPCdec;
        // Get the OTPC sensitive cell corresponding to where the track hit the sensors plane
        IDetectorElement deBarrel = calOTPCBarrel.getDetectorElement().findDetectorElement(trackPosAtSensors);
        IDetectorElement deEndcap = calOTPCEndcap.getDetectorElement().findDetectorElement(trackPosAtSensors);
        if (deBarrel == null && deEndcap != null) {
            calOTPCdec = calOTPCEndcap.getIDDecoder();
        } else if (deBarrel != null && deEndcap == null) {
            calOTPCdec = calOTPCBarrel.getIDDecoder();
        } else {
            calOTPCdec = null;
        }

        return calOTPCdec;
    }
}
